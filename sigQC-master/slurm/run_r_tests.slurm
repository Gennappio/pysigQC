#!/bin/bash
#SBATCH --job-name="sigQC_R_tests"
#SBATCH --mem=8000
#SBATCH --account=abbruzzese
#SBATCH --partition=long_gpunew
#SBATCH --cpus-per-task=4
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

set -e
echo "Job ID: $SLURM_JOB_ID  |  Node: $SLURMD_NODENAME  |  Start: $(date)"

SCRIPT_DIR="$(realpath "${SLURM_SUBMIT_DIR}")"
PROJECT_DIR="$(realpath "${SCRIPT_DIR}/..")"
ENV_NAME="sigqc_env"

# ---- Initialize conda ----
source "$(conda info --base)/etc/profile.d/conda.sh"

# ---- Create or reuse conda environment ----
if ! conda env list | grep -q "^${ENV_NAME} "; then
  echo "[INSTALL] Creating conda environment..."
  conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
    "r-base=4.3.3" "libblas=*=*openblas*" \
    r-mass r-lattice r-nnet r-class r-matrix \
    zlib libxml2 hdf5 \
    r-gridgraphics r-gplots r-fmsb r-moments r-mclust r-circlize r-biocmanager \
    r-testthat \
    bioconductor-complexheatmap bioconductor-gsva
  echo "[INSTALL] Conda environment created."
else
  echo "[INSTALL] Conda environment '${ENV_NAME}' already exists."
fi

conda activate "${ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# ---- Install biclust if needed ----
if ! Rscript -e 'library(biclust)' &>/dev/null; then
  echo "[INSTALL] Installing biclust dependencies..."
  Rscript -e 'install.packages(c("flexclust","additivityTests","tidyr","ggplot2"), repos="https://cran.r-project.org/")'
  Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/biclust/biclust_2.0.3.1.tar.gz", repos=NULL, type="source")'
else
  echo "[INSTALL] biclust already installed."
fi

# ---- Install testthat if needed ----
if ! Rscript -e 'library(testthat)' &>/dev/null; then
  Rscript -e 'install.packages("testthat", repos="https://cran.r-project.org/")'
fi

echo "[INSTALL] Done."

# ---- Run R unit tests ----
cd "${PROJECT_DIR}"

echo ""
echo "========== R UNIT TESTS (testthat) =========="
echo "Working directory: $(pwd)"
echo ""

# Step 1: Regenerate fixtures (ensures deterministic state)
echo "[FIXTURE] Regenerating test fixtures..."
Rscript tests/fixtures/fixture_generator.R
Rscript tests/fixtures/fixture_medium_generator.R
echo "[FIXTURE] Done."

# Step 2: Generate reference outputs
echo "[REFERENCE] Generating reference outputs..."
Rscript tests/fixtures/generate_reference_outputs.R
Rscript tests/fixtures/generate_medium_reference.R
echo "[REFERENCE] Done."

# Step 3: Run testthat
echo ""
echo "[TEST] Running testthat on small fixture..."
Rscript -e '
library(testthat)
library(mclust)
library(GSVA)
library(biclust)
for (f in list.files("R_refactored", pattern = "[.]R$", full.names = TRUE)) source(f)
test_dir("tests/testthat", reporter = "summary")
'

echo ""
echo "========== R TESTS COMPLETE =========="
echo "End time: $(date)"
