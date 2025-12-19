#!/usr/bin/env bash
set -euo pipefail

# ============================
# Run all 3'-end truncations of a motif through fuzznuc,
# parse outputs, and summarize hits per pattern length
# (also counting hits in first N nt of the hit start position)
#
# This version DOES NOT use FASTQ files or percentages.
# It supports .fa/.fasta and .fa.gz/.fasta.gz inputs by
# decompressing gzipped FASTA to a temp file for fuzznuc.
# ============================

# Save full command line for reporting at the top of output tables
CMDLINE="$0 $*"

# Defaults
SL_SEQ=""          # will be required via -s
MINLEN=10          # minimum length from 3' end
PMISMATCH=6
OUTDIR="fuzznuc_SL_trunc"
OVERWRITE=false
FIRST_N=30         # window for "first N nt" (hits_firstN, etc.)

usage() {
    cat <<EOF
Usage: $(basename "$0") -s SL_SEQUENCE [options] <fasta1> [fasta2 fasta3 ...]

Generate all possible patterns from the 3' end of the provided spliced leader
sequence (from length MINLEN up to full length) and run fuzznuc for each pattern
on the given FASTA file(s). For each fuzznuc output, also generate a parsed TSV,
then summarize counts across all TSVs.

Required:
  -s, --sl-seq STR      Full spliced leader sequence (5'->3')

Positional arguments:
  fasta                 One or more input FASTA files.
                        Supports: .fa/.fasta and .fa.gz/.fasta.gz

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
echo "First-N window   : $FIRST_N nt"
echo

# Optional: log all patterns used
PATTERN_LOG="${OUTDIR}/patterns_used.tsv"
echo -e "length\tpattern" > "$PATTERN_LOG"

# ----------------------------
# Handle gzipped FASTA
# ----------------------------
TMP_FASTAS=()

cleanup_tmp_fastas() {
  if (( ${#TMP_FASTAS[@]} > 0 )); then
    rm -f "${TMP_FASTAS[@]}"
  fi
}
trap cleanup_tmp_fastas EXIT

prep_fasta_for_fuzznuc() {
  local fasta="$1"
  if [[ "$fasta" == *.gz ]]; then
    local tmp
    tmp=$(mktemp "${OUTDIR}/tmp_fasta.XXXXXX.fa")
    gzip -cd "$fasta" > "$tmp"
    TMP_FASTAS+=("$tmp")
    echo "$tmp"
  else
    echo "$fasta"
  fi
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
        base_no_gz="${base%.gz}"   # strip .gz for correct stem parsing

        # sample stem from FASTA name
        stem="${base_no_gz%.fa}"
        stem="${stem%.fasta}"

        out="${OUTDIR}/${stem}_SLlen${len}.fuzznuc"
        tsv="${OUTDIR}/${stem}_SLlen${len}.tsv"

        if [[ -f "$out" && "$OVERWRITE" = false ]]; then
            echo "  Skipping fuzznuc for ${fasta} len=${len} (output exists: $out). Use -f to overwrite."
        else
            echo "  Running fuzznuc on ${fasta} -> ${out}"

            seq_in=$(prep_fasta_for_fuzznuc "$fasta")

            fuzznuc \
                -sequence "$seq_in" \
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

echo "Summaries done."
echo "Per-sample summary:    $summary_sample_len"
echo "By-length summary:     $summary_len"
echo
echo "Each summary file now starts with:"
echo "  # Command: $CMDLINE"
echo "so you always know exactly how it was generated."
