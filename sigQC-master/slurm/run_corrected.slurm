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

# ---- Install R dependencies (skipped if already installed) ----
echo "[INSTALL] Checking/installing R dependencies..."

Rscript -e '
  cran_pkgs <- c("gridGraphics","gplots","fmsb","moments","mclust","circlize","BiocManager")
  missing <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly=TRUE)]
  if (length(missing) > 0) {
    cat("Installing CRAN packages:", paste(missing, collapse=", "), "\n")
    install.packages(missing, repos="https://cran.r-project.org/")
  }

  # biclust is archived on CRAN — install from source tarball
  if (!requireNamespace("biclust", quietly=TRUE)) {
    cat("Installing biclust from CRAN archive...\n")
    tmp <- tempfile(fileext=".tar.gz")
    download.file("https://cran.r-project.org/src/contrib/Archive/biclust/biclust_2.0.3.1.tar.gz", tmp)
    install.packages(tmp, repos=NULL, type="source")
  }

  # Bioconductor packages
  bioc_pkgs <- c("ComplexHeatmap","GSVA")
  missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly=TRUE)]
  if (length(missing_bioc) > 0) {
    cat("Installing Bioconductor packages:", paste(missing_bioc, collapse=", "), "\n")
    BiocManager::install(missing_bioc, ask=FALSE)
  }

  cat("All dependencies OK.\n")
'

# ---- Install sigQC from source (skipped if already installed) ----
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
