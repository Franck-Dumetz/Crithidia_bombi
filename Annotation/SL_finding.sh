#!/usr/bin/env bash
set -euo pipefail

# ============================
# Run all 3'-end truncations of a motif through fuzznuc,
# parse outputs, and summarize hits per pattern length
# with percentages based on FASTQ read counts
# (also counting hits in first N nt of the read)
# ============================

# Save full command line for reporting at the top of output tables
CMDLINE="$0 $*"

# Defaults
SL_SEQ=""          # will be required via -s
MINLEN=10          # minimum length from 3' end
PMISMATCH=6
OUTDIR="fuzznuc_SL_trunc"
OVERWRITE=false
FASTQ_DIR=""       # will be required for percentages
FIRST_N=30         # window for "first N nt" (hits_firstN, etc.)

usage() {
    cat <<EOF
Usage: $(basename "$0") -s SL_SEQUENCE -F FASTQ_DIR [options] <fasta1> [fasta2 fasta3 ...]

Generate all possible patterns from the 3' end of the provided spliced leader
sequence (from length MINLEN up to full length) and run fuzznuc for each pattern
on the given FASTA file(s). For each fuzznuc output, also generate a parsed TSV,
then summarize counts across all TSVs and compute percentages based on the number
of reads in the corresponding FASTQ file.

Required:
  -s, --sl-seq STR      Full spliced leader sequence (5'->3')
  -F, --fastq-dir DIR   Directory containing FASTQ files for the samples.
                        FASTQ filenames must share the same stem as the FASTA,
                        e.g. sample.fa -> sample.fastq.gz / sample.fastq / sample.fq.gz / sample.fq

Positional arguments:
  fasta                 One or more input FASTA files.

Options:
  -l, --minlen INT      Minimum pattern length (from 3' end; default: ${MINLEN})
  -m, --mismatch INT    Number of allowed mismatches (pmismatch; default: ${PMISMATCH})
  -o, --outdir DIR      Output directory for all fuzznuc results (default: ${OUTDIR})
  -n, --first-n INT     Window from 5' end to count hits (default: ${FIRST_N})
  -f, --force           Overwrite existing output files
  -h, --help            Show this help and exit
EOF
}

# ---- Parse arguments ----
ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--sl-seq)
            SL_SEQ="$2"
            shift 2
            ;;
        -l|--minlen)
            MINLEN="$2"
            shift 2
            ;;
        -m|--mismatch)
            PMISMATCH="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -n|--first-n)
            FIRST_N="$2"
            shift 2
            ;;
        -f|--force)
            OVERWRITE=true
            shift
            ;;
        -F|--fastq-dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

# Add any remaining positional args, if present
if [[ $# -gt 0 ]]; then
    ARGS+=("$@")
fi

# ---- Check inputs ----
if [[ -z "$SL_SEQ" ]]; then
    echo "ERROR: You must provide the spliced leader sequence with -s/--sl-seq." >&2
    usage
    exit 1
fi

if [[ -z "$FASTQ_DIR" ]]; then
    echo "ERROR: You must provide the FASTQ directory with -F/--fastq-dir." >&2
    usage
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR" >&2
    exit 1
fi

if [[ ${#ARGS[@]} -eq 0 ]]; then
    echo "ERROR: You must provide at least one FASTA file." >&2
    usage
    exit 1
fi

# Normalize SL sequence (uppercase, remove spaces)
SL_SEQ=$(echo "$SL_SEQ" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')
SL_LEN=${#SL_SEQ}

if (( MINLEN < 1 || MINLEN > SL_LEN )); then
    echo "ERROR: Invalid --minlen (${MINLEN}); must be between 1 and SL length (${SL_LEN})." >&2
    exit 1
fi

# Check FIRST_N is a positive integer
if ! [[ "$FIRST_N" =~ ^[0-9]+$ ]] || (( FIRST_N <= 0 )); then
    echo "ERROR: --first-n / -n must be a positive integer (got: $FIRST_N)" >&2
    exit 1
fi

# Check fuzznuc is available
if ! command -v fuzznuc >/dev/null 2>&1; then
    echo "ERROR: fuzznuc not found in PATH. Load the EMBOSS module or install it." >&2
    exit 1
fi

mkdir -p "$OUTDIR"

echo "Full SL sequence : $SL_SEQ"
echo "SL length        : $SL_LEN"
echo "Min pattern len  : $MINLEN"
echo "Mismatches       : $PMISMATCH"
echo "Output dir       : $OUTDIR"
echo "FASTA files      : ${ARGS[*]}"
echo "FASTQ dir        : $FASTQ_DIR"
echo "First-N window   : $FIRST_N nt"
echo

# Optional: log all patterns used
PATTERN_LOG="${OUTDIR}/patterns_used.tsv"
echo -e "length\tpattern" > "$PATTERN_LOG"

# File for total reads per sample
READS_TSV="${OUTDIR}/reads_per_sample.tsv"
echo -e "sample\treads" > "$READS_TSV"

# Track reads per sample in-memory too
declare -A READS_PER_SAMPLE

# Helper: find FASTQ file for a given sample stem
find_fastq_for_sample() {
    local stem="$1"
    local fq=""
    local d="$FASTQ_DIR"

    local candidates=(
        "$d/${stem}.fastq.gz"
        "$d/${stem}.fastq"
        "$d/${stem}.fq.gz"
        "$d/${stem}.fq"
    )

    for c in "${candidates[@]}"; do
        if [[ -f "$c" ]]; then
            fq="$c"
            break
        fi
    done

    if [[ -z "$fq" ]]; then
        echo "ERROR: No FASTQ found for sample '$stem' in $FASTQ_DIR" >&2
        exit 1
    fi

    echo "$fq"
}

# Helper: count reads in FASTQ (lines/4)
count_reads_in_fastq() {
    local fq="$1"
    local lines=0
    if [[ "$fq" == *.gz ]]; then
        lines=$(zcat "$fq" | wc -l)
    else
        lines=$(wc -l < "$fq")
    fi
    echo $(( lines / 4 ))
}

# ---- Main: loop over pattern lengths and FASTA files ----
for (( len = MINLEN; len <= SL_LEN; len++ )); do
    # Take the last 'len' nt (3' end) of SL_SEQ
    start=$(( SL_LEN - len ))
    pattern="${SL_SEQ:start:len}"

    echo -e "${len}\t${pattern}" >> "$PATTERN_LOG"
    echo "=== Pattern length ${len} nt (3'-end) ==="
    echo "Pattern: ${pattern}"

    for fasta in "${ARGS[@]}"; do
        if [[ ! -f "$fasta" ]]; then
            echo "  WARNING: FASTA not found, skipping: $fasta" >&2
            continue
        fi

        base=$(basename "$fasta")
        # e.g. monocistron.fa -> monocistron_SLlen10.fuzznuc / .tsv
        stem="${base%.fa}"
        stem="${stem%.fasta}"

        # Ensure we have read count for this sample
        if [[ -z "${READS_PER_SAMPLE[$stem]:-}" ]]; then
            fq=$(find_fastq_for_sample "$stem")
            echo "  Counting reads in FASTQ for sample '$stem': $fq"
            reads=$(count_reads_in_fastq "$fq")
            READS_PER_SAMPLE[$stem]=$reads
            echo -e "${stem}\t${reads}" >> "$READS_TSV"
            echo "  -> $reads reads"
        fi

        out="${OUTDIR}/${stem}_SLlen${len}.fuzznuc"
        tsv="${OUTDIR}/${stem}_SLlen${len}.tsv"

        if [[ -f "$out" && "$OVERWRITE" = false ]]; then
            echo "  Skipping fuzznuc for ${fasta} len=${len} (output exists: $out). Use -f to overwrite."
        else
            echo "  Running fuzznuc on ${fasta} -> ${out}"
            fuzznuc \
                -sequence "$fasta" \
                -pattern "$pattern" \
                -pmismatch "$PMISMATCH" \
                -outfile "$out"
        fi

        # ---- Parse fuzznuc output into compact TSV (1 line per hit) ----
        if [[ -f "$out" ]]; then
            echo "  Parsing ${out} -> ${tsv}"
            grep 'Sequence' --no-group-separator -A1 "$out" \
                | grep -v 'HitCount' \
                | grep -v 'Start' \
                | awk '{printf "%s%s",$0,NR%2?"\t":RS}' \
                > "$tsv"
        else
            echo "  WARNING: Expected fuzznuc output not found for len=${len}: $out" >&2
        fi
    done

    echo
done

echo "All patterns done."
echo "Pattern list written to: $PATTERN_LOG"
echo "Reads per sample:        $READS_TSV"
echo "fuzznuc outputs and TSVs in: $OUTDIR/"
echo

# ============================
# Summarize TSVs
# ============================

echo "Summarizing TSV files..."

summary_sample_len="${OUTDIR}/summary_by_sample_and_length.tsv"
summary_len="${OUTDIR}/summary_by_length.tsv"

echo "Writing per-sample summary to: $summary_sample_len"
echo "Writing aggregated summary to:  $summary_len"

# Header for per-sample summary
echo "# Command: $CMDLINE" > "$summary_sample_len"
echo -e "sample\tlength\thits\thits_first${FIRST_N}" >> "$summary_sample_len"

shopt -s nullglob
for f in "$OUTDIR"/*_SLlen*.tsv; do
    fname=$(basename "$f")

    # sample = part before "_SLlen"
    sample="${fname%%_SLlen*}"

    # length = number between SLlen and .tsv
    len_part="${fname#*_SLlen}"  # e.g. "10.tsv"
    len="${len_part%.tsv}"       # e.g. "10"

    # total hits = number of lines
    hits=$(wc -l < "$f")

    # hits with start <= FIRST_N
    hits_firstN=$(awk -v maxpos="$FIRST_N" '
        {
            # Split into "Sequence..." and "start end strand pattern..." parts
            n = split($0, a, "\t")
            if (n < 2) next
            split(a[2], b, " ")
            start = b[1] + 0
            if (start > 0 && start <= maxpos) c++
        }
        END { print c+0 }
    ' "$f")

    echo -e "${sample}\t${len}\t${hits}\t${hits_firstN}" >> "$summary_sample_len"
done
shopt -u nullglob

# Aggregate by length (sum across samples)
echo "# Command: $CMDLINE" > "$summary_len"
echo -e "length\thits\thits_first${FIRST_N}" >> "$summary_len"
awk '
    FNR <= 2 { next }   # skip command + header lines in summary_sample_len
    {
        len  = $2
        h    = $3
        hN   = $4
        total[len]  += h
        totalN[len] += hN
    }
    END {
        for (l in total) {
            printf "%s\t%s\t%s\n", l, total[l], totalN[l]
        }
    }
' "$summary_sample_len" | sort -n >> "$summary_len"

# ============================
# Add percentage columns based on FASTQ read counts
# ============================

echo "Adding percentage columns based on FASTQ read counts..."

# Per-sample + length:
#   percent_of_sample_reads           = hits / sample_reads
#   percent_of_sample_reads_firstN    = hits_firstN / sample_reads
tmp="${summary_sample_len}.tmp"
{
    # We want the command line to remain first in the *new* file
    echo "# Command: $CMDLINE"
    awk -v N="$FIRST_N" '
        NR==FNR && FNR>1 {
            reads[$1]=$2      # from reads_per_sample.tsv: sample -> total reads
            next
        }
        NR==FNR { next }      # header of reads_per_sample.tsv

        FNR==1 { next }       # old "# Command" line in summary_by_sample_and_length
        FNR==2 {
            print $0 "\tpercent_of_sample_reads\tpercent_of_sample_reads_first" N
            next
        }
        {
            samp = $1
            r    = reads[samp]
            h    = $3
            hN   = $4
            pct  = (r > 0 ? (h  * 100.0 / r) : 0)
            pctN = (r > 0 ? (hN * 100.0 / r) : 0)
            printf "%s\t%.6f\t%.6f\n", $0, pct, pctN
        }
    ' "$READS_TSV" "$summary_sample_len"
} > "$tmp"
mv "$tmp" "$summary_sample_len"

# By length:
#   percent_of_all_reads           = hits / total_reads_all
#   percent_of_all_reads_firstN    = hits_firstN / total_reads_all
total_reads_all=$(awk 'NR>1 {s+=$2} END {print s+0}' "$READS_TSV")

tmp="${summary_len}.tmp"
{
    echo "# Command: $CMDLINE"
    awk -v tot="$total_reads_all" -v N="$FIRST_N" '
        FNR==1 { next }   # old "# Command"
        FNR==2 {
            print $0 "\tpercent_of_all_reads\tpercent_of_all_reads_first" N
            next
        }
        {
            h  = $2
            hN = $3
            pct  = (tot > 0 ? (h  * 100.0 / tot) : 0)
            pctN = (tot > 0 ? (hN * 100.0 / tot) : 0)
            printf "%s\t%.6f\t%.6f\n", $0, pct, pctN
        }
    ' "$summary_len"
} > "$tmp"
mv "$tmp" "$summary_len"

echo "Summaries done."
echo "Per-sample summary:    $summary_sample_len"
echo "By-length summary:     $summary_len"
echo
echo "Each summary file now starts with:"
echo "  # Command: $CMDLINE"
echo "so you always know exactly how it was generated."
