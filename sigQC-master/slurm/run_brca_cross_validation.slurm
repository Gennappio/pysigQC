#!/bin/bash
#SBATCH --job-name="sigQC_brca_xval"
#SBATCH --mem=16000
#SBATCH --account=abbruzzese
#SBATCH --partition=long_gpunew
#SBATCH --cpus-per-task=6
#SBATCH --output=out/%x_%j.out
#SBATCH --error=err/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=gennaro.abbruzzese@unibocconi.it

set -e
echo "Job ID: $SLURM_JOB_ID  |  Node: $SLURMD_NODENAME  |  Start: $(date)"

SCRIPT_DIR="$(realpath "${SLURM_SUBMIT_DIR}")"
PROJECT_DIR="$(realpath "${SCRIPT_DIR}/../..")"
R_ENV_NAME="sigqc_env"
PY_ENV_NAME="sigqc_py_env"

INPUT_DIR="${PROJECT_DIR}/sigQC-master/input_files"
BRCA_XVAL_DIR="${PROJECT_DIR}/brca_cross_validation"

if [ ! -d "${INPUT_DIR}" ]; then
  echo "[ERROR] BRCA input_files not found at ${INPUT_DIR}."
  echo "        Place the original paper data (expression matrices + signatures) there."
  exit 1
fi

mkdir -p "${BRCA_XVAL_DIR}/r_output" "${BRCA_XVAL_DIR}/py_output"

# ============================================================================
# STEP 1: Run R refactored pipeline on BRCA data
# ============================================================================
echo ""
echo "========== STEP 1: R REFACTORED PIPELINE ON BRCA DATA =========="

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${R_ENV_NAME}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

cd "${INPUT_DIR}"

cat > "${BRCA_XVAL_DIR}/run_r_brca.R" << 'RSCRIPT'
library(mclust)
library(GSVA)
library(biclust)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- args[1]
out_dir <- args[2]

# Source refactored R code
for (f in list.files(file.path(project_dir, "sigQC-master", "R_refactored"),
                     pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

# Load BRCA data (same as the original paper vignette)
# Adjust these file names if your input_files/ layout differs
input_dir <- file.path(project_dir, "sigQC-master", "input_files")

# Discover expression matrices and signatures
# Expects: expression_*.csv and signatures.csv (or .rds)
expr_files <- list.files(input_dir, pattern = "^exprs.*\\.(csv|rds)$", full.names = TRUE)
sig_files  <- list.files(input_dir, pattern = "^sig.*\\.(csv|rds)$", full.names = TRUE)

if (length(expr_files) == 0 || length(sig_files) == 0) {
  cat("[ERROR] Could not find expression/signature files in", input_dir, "\n")
  cat("Expression files found:", paste(basename(expr_files), collapse=", "), "\n")
  cat("Signature files found:", paste(basename(sig_files), collapse=", "), "\n")
  quit(status = 1)
}

# Load expression matrices
mRNA_expr_matrix <- list()
for (f in expr_files) {
  ds_name <- tools::file_path_sans_ext(basename(f))
  if (grepl("\\.rds$", f)) {
    mRNA_expr_matrix[[ds_name]] <- readRDS(f)
  } else {
    mRNA_expr_matrix[[ds_name]] <- as.matrix(read.csv(f, row.names = 1))
  }
}
names_datasets <- names(mRNA_expr_matrix)

# Load signatures
sig_file <- sig_files[1]
if (grepl("\\.rds$", sig_file)) {
  gene_sigs_list <- readRDS(sig_file)
} else {
  sig_df <- read.csv(sig_file)
  gene_sigs_list <- list()
  for (sn in unique(sig_df$signature)) {
    genes <- sig_df$gene[sig_df$signature == sn]
    gene_sigs_list[[sn]] <- as.matrix(genes)
    colnames(gene_sigs_list[[sn]]) <- sn
  }
}
names_sigs <- names(gene_sigs_list)

cat("Datasets:", paste(names_datasets, collapse=", "), "\n")
cat("Signatures:", paste(names_sigs, collapse=", "), "\n")

# Run all compute modules
cat("Running compute_var...\n")
var_r <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)

cat("Running compute_expr...\n")
expr_r <- compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)

cat("Running compute_compactness...\n")
compact_r <- compute_compactness(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)

cat("Running compute_stan...\n")
stan_r <- compute_stan(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)

cat("Running compute_metrics...\n")
metrics_r <- compute_metrics(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)

# Save per-module radar values
for (sig in names_sigs) {
  for (ds in names_datasets) {
    for (mod_name in c("var", "expr", "compact", "stan", "metrics")) {
      mod <- switch(mod_name,
                    var = var_r, expr = expr_r, compact = compact_r,
                    stan = stan_r, metrics = metrics_r)
      rv <- mod$radar_values[[sig]][[ds]]
      write.csv(t(as.data.frame(rv)),
                file.path(out_dir, paste0(mod_name, "_radar_", sig, "_", ds, ".csv")),
                row.names = FALSE)
    }
  }
}

# Assemble full radar
radar_values <- list()
for (sig in names_sigs) {
  radar_values[[sig]] <- list()
  for (ds in names_datasets) {
    radar_values[[sig]][[ds]] <- c(
      var_r$radar_values[[sig]][[ds]],
      expr_r$radar_values[[sig]][[ds]],
      compact_r$radar_values[[sig]][[ds]],
      metrics_r$radar_values[[sig]][[ds]],
      stan_r$radar_values[[sig]][[ds]]
    )
  }
}
radar_result <- compute_radar(radar_values, names_sigs, names_datasets)
write.csv(radar_result$output_table,
          file.path(out_dir, "radar_output_table.csv"))

cat("R BRCA outputs saved to", out_dir, "\n")
RSCRIPT

Rscript "${BRCA_XVAL_DIR}/run_r_brca.R" "${PROJECT_DIR}" "${BRCA_XVAL_DIR}/r_output"

# ============================================================================
# STEP 2: Run Python pipeline on BRCA data
# ============================================================================
echo ""
echo "========== STEP 2: PYTHON PIPELINE ON BRCA DATA =========="

conda activate "${PY_ENV_NAME}"
cd "${PROJECT_DIR}"

cat > "${BRCA_XVAL_DIR}/run_py_brca.py" << 'PYSCRIPT'
"""Run Python pipeline on BRCA data and save outputs for cross-validation."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

project_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
r_out_dir = Path(sys.argv[3])

sys.path.insert(0, str(project_dir))

from pysigqc.eval_var import compute_var
from pysigqc.eval_expr import compute_expr
from pysigqc.eval_compactness import compute_compactness
from pysigqc.eval_stan import compute_stan
from pysigqc.compare_metrics import compute_metrics
from pysigqc.radar_chart import compute_radar

input_dir = project_dir / "sigQC-master" / "input_files"

# Load expression matrices
import glob
expr_files = sorted(input_dir.glob("exprs*.csv")) + sorted(input_dir.glob("exprs*.rds"))
datasets = {}
for f in expr_files:
    ds_name = f.stem
    if f.suffix == ".csv":
        datasets[ds_name] = pd.read_csv(f, index_col=0)
    else:
        print(f"[WARN] Cannot read .rds in Python: {f}. Convert to CSV first.")
        continue
names_datasets = list(datasets.keys())

# Load signatures
sig_files = sorted(input_dir.glob("sig*.csv"))
if not sig_files:
    print("[ERROR] No signature CSV files found")
    sys.exit(1)
sig_df = pd.read_csv(sig_files[0])
signatures = {}
for sig_name, group in sig_df.groupby("signature"):
    signatures[sig_name] = group["gene"].tolist()
names_sigs = list(signatures.keys())

print(f"Datasets: {names_datasets}")
print(f"Signatures: {names_sigs}")

# Run all compute modules
print("Running compute_var...")
var_r = compute_var(signatures, names_sigs, datasets, names_datasets)

print("Running compute_expr...")
expr_r = compute_expr(signatures, names_sigs, datasets, names_datasets)

print("Running compute_compactness...")
compact_r = compute_compactness(signatures, names_sigs, datasets, names_datasets)

print("Running compute_stan...")
stan_r = compute_stan(signatures, names_sigs, datasets, names_datasets)

print("Running compute_metrics...")
metrics_r = compute_metrics(signatures, names_sigs, datasets, names_datasets)

# Assemble full radar
radar_values = {}
for sig in names_sigs:
    radar_values[sig] = {}
    for ds in names_datasets:
        vals = {}
        vals.update(var_r["radar_values"][sig][ds])
        vals.update(expr_r["radar_values"][sig][ds])
        vals.update(compact_r["radar_values"][sig][ds])
        vals.update(metrics_r["radar_values"][sig][ds])
        vals.update(stan_r["radar_values"][sig][ds])
        radar_values[sig][ds] = vals

radar_result = compute_radar(radar_values, names_sigs, names_datasets)
radar_result["output_table"].to_csv(out_dir / "radar_output_table.csv")

print(f"Python BRCA outputs saved to {out_dir}")

# ============================================================================
# STEP 3: Compare R vs Python outputs
# ============================================================================
print("\n========== CROSS-VALIDATION: R vs Python ==========")

r_table = pd.read_csv(r_out_dir / "radar_output_table.csv", index_col=0)
py_table = radar_result["output_table"]

# Align indices
common_idx = r_table.index.intersection(py_table.index)
common_cols = r_table.columns.intersection(py_table.columns)

print(f"Comparing {len(common_idx)} rows x {len(common_cols)} columns")

mismatches = []
for col in common_cols:
    for idx in common_idx:
        r_val = r_table.loc[idx, col]
        py_val = py_table.loc[idx, col]
        if np.isnan(r_val) and np.isnan(py_val):
            continue
        if np.isnan(r_val) or np.isnan(py_val):
            mismatches.append((idx, col, r_val, py_val, "NaN mismatch"))
            continue
        # PCA sign ambiguity
        if "pca1" in col:
            rdiff = abs(abs(py_val) - abs(r_val))
            rtol_check = abs(r_val) * 0.05
        else:
            rdiff = abs(py_val - r_val)
            rtol_check = abs(r_val) * 0.01
        if rdiff > max(rtol_check, 1e-8):
            mismatches.append((idx, col, r_val, py_val, f"diff={rdiff:.6f}"))

if mismatches:
    print(f"\n[WARN] {len(mismatches)} mismatches found:")
    mismatch_df = pd.DataFrame(mismatches, columns=["row", "col", "R", "Python", "note"])
    print(mismatch_df.to_string(index=False))
    mismatch_df.to_csv(out_dir / "brca_mismatches.csv", index=False)
else:
    print("\n[OK] All values match within tolerance.")

print("\nCross-validation complete.")
PYSCRIPT

python "${BRCA_XVAL_DIR}/run_py_brca.py" "${PROJECT_DIR}" "${BRCA_XVAL_DIR}/py_output" "${BRCA_XVAL_DIR}/r_output"

echo ""
echo "========== BRCA CROSS-VALIDATION COMPLETE =========="
echo "R outputs:      ${BRCA_XVAL_DIR}/r_output/"
echo "Python outputs: ${BRCA_XVAL_DIR}/py_output/"
echo "End time: $(date)"
