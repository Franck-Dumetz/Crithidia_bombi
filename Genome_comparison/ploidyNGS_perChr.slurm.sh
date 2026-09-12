#!/usr/bin/env bash
#SBATCH --job-name=ploidyNGS_chr
#SBATCH --output=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/ploidy/ploidyNGS_perChr/logs/%x.%A_%a.out
#SBATCH --error=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/ploidy/ploidyNGS_perChr/logs/%x.%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --account=dserre-lab
#SBATCH --time=06:00:00
#SBATCH --array=1-34
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fdumetz@som.umaryland.edu

set -euo pipefail

# -------------------------
# Paths
# -------------------------
BAM_DIR="/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/ploidy/ploidyNGS_perChr"
OUT_DIR="/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/ploidy/ploidyNGS_perChr"

PLOIDYNGS_DIR="/local/projects-t3/SerreDLab-3/fdumetz/programs/ploidyNGS"
PLOIDYNGS_PY="${PLOIDYNGS_DIR}/ploidyNGS.py"

# -------------------------
# ploidyNGS parameters
# -------------------------
MINCOV=20       # -u
MAXDEPTH=200    # -d
MAXALLELE=0.95  # -m

mkdir -p "$OUT_DIR"/{results,logs}

# -------------------------
# Pick BAM for this array task
# -------------------------
mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)
N="${#BAMS[@]}"

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [[ "$IDX" -lt 0 || "$IDX" -ge "$N" ]]; then
  echo "ERROR: Array index out of range. Task ${SLURM_ARRAY_TASK_ID} but only ${N} BAMs."
  exit 1
fi

BAM="${BAMS[$IDX]}"
BASE="$(basename "$BAM" .bam)"
PREFIX="$OUT_DIR/results/${BASE}"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}/${N}"
echo "BAM:    $BAM"
echo "PREFIX: $PREFIX"

# -------------------------
# Activate conda env + run ploidyNGS
# -------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ploidyngs

cd "$PLOIDYNGS_DIR"
python "$PLOIDYNGS_PY" -o "$PREFIX" -b "$BAM" -u "$MINCOV" -d "$MAXDEPTH" -m "$MAXALLELE"

conda deactivate
echo "[$(date)] Done: $BASE"
