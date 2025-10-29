#!/bin/bash
#SBATCH --job-name=Cb_aligned
#SBATCH --output=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/DGE/CbHisat.out
#SBATCH --error=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/DGE/CbHisat.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --array=1-6

set -euo pipefail

# -------- CONFIG --------
HISAT2_INDEX="/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/DGE/Cb_WHA1_index/Cbombi_WHA1_ref"
FASTQ_ROOT="/local/projects-t3/EDS10"                 # where we search
ALIGN_OUTPUT_DIR="/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/DGE/bam_files"
SUM_DIR="${ALIGN_OUTPUT_DIR}/sum_dir"
THREADS="${SLURM_CPUS_PER_TASK:-12}"
PICARD="/usr/local/packages/picard-2.9.4/picard.jar"

mkdir -p "$ALIGN_OUTPUT_DIR" "$SUM_DIR"

# -------- DISCOVER PAIRED TRIMMED FASTQs --------
declare -a SAMPLE_TABLE
while IFS= read -r -d '' R1; do
  # Expect matching R2 in same dir
  R2="${R1/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}"
  if [[ -f "$R2" ]]; then
    base="$(basename "$R1")"
    sample="${base%_R1_trimmed.fastq.gz}"       # use long prefix as sample id
    SAMPLE_TABLE+=("${sample}"$'\t'"${R1}"$'\t'"${R2}")
  else
    echo "WARNING: Missing R2 for $R1" >&2
  fi
done < <(find "$FASTQ_ROOT" -type f -name "*_R1_trimmed.fastq.gz" -path "*/Cbombi_WHA1*/ILLUMINA_DATA/*" -print0)

NUM_SAMPLES="${#SAMPLE_TABLE[@]}"
if (( NUM_SAMPLES == 0 )); then
  echo "ERROR: No trimmed pairs found under $FASTQ_ROOT/*/Cbombi_WHA1*/ILLUMINA_DATA/." >&2
  exit 2
fi

TASK_INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))
if (( TASK_INDEX < 0 || TASK_INDEX >= NUM_SAMPLES )); then
  echo "No work for array task ${SLURM_ARRAY_TASK_ID} (found ${NUM_SAMPLES} pairs)."
  exit 0
fi

IFS=$'\t' read -r SAMPLE FILE1 FILE2 <<< "${SAMPLE_TABLE[$TASK_INDEX]}"

echo "Pairs found: ${NUM_SAMPLES} | This task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: $SAMPLE"
echo "R1: $FILE1"
echo "R2: $FILE2"

# -------- ALIGN --------
/usr/local/packages/hisat2-2.2.1/hisat2 \
  -x "$HISAT2_INDEX" \
  -1 "$FILE1" -2 "$FILE2" \
  -S "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam" \
  --max-intronlen 5000 \
  -p "$THREADS" \
  --summary-file "$SUM_DIR/${SAMPLE}_summary.txt" \
  --new-summary

# Convert/filter/sort
samtools view -@ "$THREADS" -bhF 2308 "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam" > "$ALIGN_OUTPUT_DIR/${SAMPLE}.bam"
samtools sort -@ "$THREADS" "$ALIGN_OUTPUT_DIR/${SAMPLE}.bam" -o "$ALIGN_OUTPUT_DIR/${SAMPLE}_sorted.bam"

# -------- DEDUP --------
java -Xmx44g -jar "$PICARD" MarkDuplicates \
  INPUT="$ALIGN_OUTPUT_DIR/${SAMPLE}_sorted.bam" \
  OUTPUT="$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup.bam" \
  METRICS_FILE="$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup_metrics.txt" \
  REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
  TMP_DIR="${TMPDIR:-/tmp}"

[[ -s "$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup.bam" ]] || { echo "ERROR: Picard failed for $SAMPLE"; exit 3; }
samtools index "$ALIGN_OUTPUT_DIR/${SAMPLE}_dedup.bam"

# Cleanup
rm -f "$ALIGN_OUTPUT_DIR/${SAMPLE}.sam" "$ALIGN_OUTPUT_DIR/${SAMPLE}.bam" "$ALIGN_OUTPUT_DIR/${SAMPLE}_sorted.bam"

echo "Done: $SAMPLE → ${ALIGN_OUTPUT_DIR}/${SAMPLE}_dedup.bam"