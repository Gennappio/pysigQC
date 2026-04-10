#!/bin/bash
#SBATCH --job-name="sigQC_corrected"
#SBATCH --mem=16000
#SBATCH --account=abbruzzese
#SBATCH --partition=long_gpunew
#SBATCH --cpus-per-task=6
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

echo "start"
echo "Running corrected sigQC (8 bug fixes applied) on BRCA paper example data"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"

SCRIPT_DIR="$(dirname "$0")"
ENV_NAME="sigqc_env"

# ---- Initialize conda ----
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---- Create conda environment if it doesn't exist ----
if ! conda env list | grep -q "^${ENV_NAME} "; then
  echo "[INSTALL] Creating conda environment '${ENV_NAME}'..."
  conda create -y -n "${ENV_NAME}" \
    -c conda-forge -c bioconda \
    r-base \
    r-gridgraphics \
    r-biclust \
    r-gplots \
    r-fmsb \
    r-moments \
    r-mclust \
    r-circlize \
    r-biocmanager \
    bioconductor-complexheatmap \
    bioconductor-gsva
  echo "[INSTALL] Conda environment created."
else
  echo "[INSTALL] Conda environment '${ENV_NAME}' already exists, skipping creation."
fi

conda activate "${ENV_NAME}"

# ---- Install sigQC from source into the conda R ----
# Note: run_corrected_with_timing.R sources R_corrected/ directly, so the installed
# package is only needed for its namespace; the corrected functions override it at runtime.
if ! Rscript -e 'library(sigQC)' &>/dev/null; then
  echo "[INSTALL] Installing sigQC from source..."
  R CMD INSTALL "${SCRIPT_DIR}/.."
else
  echo "[INSTALL] sigQC already installed, skipping."
fi

echo "[INSTALL] Done."

# ---- Run analysis ----
cd "${SCRIPT_DIR}/../input_files"

Rscript run_corrected_with_timing.R

echo "End time: $(date)"
echo "done"
