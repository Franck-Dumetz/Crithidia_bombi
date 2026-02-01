#!/usr/bin/env bash
#SBATCH --job-name=ploidyNGS_chr
#SBATCH --output=/logs/%x.%A_%a.out
#SBATCH --error=/logs/%x.%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --array=1-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fdumetz@som.umaryland.edu

set -euo pipefail

# -------------------------
# Paths
# -------------------------
BAM_DIR="/ploidyNGS_perChr"
OUT_DIR="/ploidyNGS_perChr"

PLOIDYNGS_DIR="/ploidyNGS"
VENV_ACTIVATE="/.venv/bin/activate"
PLOIDYNGS_PY="/ploidyNGS.py"   # <-- FIXED (was pointing to the .sh)

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# -------------------------
# ploidyNGS parameters
# -------------------------
MINCOV=20       # -u
MAXDEPTH=200    # -d
MAXALLELE=0.95  # -m
DO_GUESS=0      # set to 1 to add -g

# Ensure these exist BEFORE submission too (Slurm opens logs early)
mkdir -p "$OUT_DIR"/{results,logs}

# -------------------------
# Pick BAM for this array task
# -------------------------
mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)
N="${#BAMS[@]}"

if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No BAMs found in $BAM_DIR"
  exit 1
fi

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
# Activate env + run ploidyNGS
# -------------------------
if [[ ! -f "$VENV_ACTIVATE" ]]; then
  echo "ERROR: Cannot find venv activate script: $VENV_ACTIVATE"
  exit 1
fi
if [[ ! -f "$PLOIDYNGS_PY" ]]; then
  echo "ERROR: Cannot find ploidyNGS.py at: $PLOIDYNGS_PY"
  exit 1
fi

source "$VENV_ACTIVATE"
cd "$PLOIDYNGS_DIR"

CMD=(python "$PLOIDYNGS_PY" -o "$PREFIX" -b "$BAM" -u "$MINCOV" -d "$MAXDEPTH" -m "$MAXALLELE")
if [[ "$DO_GUESS" -eq 1 ]]; then
  CMD+=(-g)
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

deactivate
echo "[$(date)] Done: $BASE"
