#!/bin/bash
#SBATCH --job-name="sigQC_original"
#SBATCH --mem=16000
#SBATCH --account=abbruzzese
#SBATCH --partition=long_gpunew
#SBATCH --cpus-per-task=6
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

set -e
echo "start"
echo "Job ID: $SLURM_JOB_ID  |  Node: $SLURMD_NODENAME  |  Start: $(date)"

SCRIPT_DIR="$(realpath "${SLURM_SUBMIT_DIR}")"
ENV_NAME="sigqc_env"

# ---- Initialize conda ----
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---- Create conda environment if it doesn't exist ----
if ! conda env list | grep -q "^${ENV_NAME} "; then
  echo "[INSTALL] Creating conda environment..."
  conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
    "r-base=4.3.3" "libblas=*=*openblas*" \
    r-mass r-lattice r-nnet r-class r-matrix \
    zlib libxml2 hdf5 \
    r-gridgraphics r-gplots r-fmsb r-moments r-mclust r-circlize r-biocmanager \
    bioconductor-complexheatmap bioconductor-gsva
  echo "[INSTALL] Conda environment created."
else
  echo "[INSTALL] Conda environment '${ENV_NAME}' already exists, skipping."
fi

conda activate "${ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# ---- Install biclust and its CRAN dependencies (skipped if present) ----
if ! Rscript -e 'library(biclust)' &>/dev/null; then
  echo "[INSTALL] Installing biclust dependencies..."
  Rscript -e 'install.packages(c("flexclust","additivityTests","tidyr","ggplot2"), repos="https://cran.r-project.org/")'
  echo "[INSTALL] Installing biclust from CRAN archive..."
  Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/biclust/biclust_2.0.3.1.tar.gz", repos=NULL, type="source")'
else
  echo "[INSTALL] biclust already installed, skipping."
fi

# ---- Install sigQC from source (skipped if present) ----
if ! Rscript -e 'library(sigQC)' &>/dev/null; then
  echo "[INSTALL] Installing sigQC..."
  R CMD INSTALL "${SCRIPT_DIR}/.."
else
  echo "[INSTALL] sigQC already installed, skipping."
fi

echo "[INSTALL] Done."

# ---- Run analysis ----
INPUT_DIR="${SCRIPT_DIR}/../input_files"
echo "[RUN] Changing to: $(realpath ${INPUT_DIR})"
cd "${INPUT_DIR}"

Rscript run_original_with_timing.R

echo "End time: $(date)"
echo "done"
