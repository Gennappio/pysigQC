#!/usr/bin/env Rscript
#
# Equivalent of:
#   python main.py --datasets "sigQC-master/tests/fixtures/fixture_dataset_*.csv" \
#                  --signatures sigQC-master/tests/fixtures/fixture_signatures.csv \
#                  --out-dir results_r/
#
# Usage:
#   Rscript run_r_sigqc.R

# --- 1. Load datasets ---
dataset_files <- Sys.glob("sigQC-master/tests/fixtures/fixture_medium_dataset_*.csv")
cat("Loading", length(dataset_files), "dataset(s):\n")

mRNA_expr_matrix <- list()
for (f in dataset_files) {
  # Dataset name = filename without extension (e.g. "fixture_dataset_A")
  name <- tools::file_path_sans_ext(basename(f))
  cat("  ", f, " -> '", name, "'\n", sep = "")
  mRNA_expr_matrix[[name]] <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
}
names_datasets <- names(mRNA_expr_matrix)

# --- 2. Load signatures (long-format CSV: signature, gene) ---
sigs_df <- read.csv("sigQC-master/tests/fixtures/fixture_medium_signatures.csv",
                     stringsAsFactors = FALSE)
gene_sigs_list <- list()
for (sig_name in unique(sigs_df$signature)) {
  genes <- sigs_df$gene[sigs_df$signature == sig_name]
  gene_sigs_list[[sig_name]] <- as.matrix(genes)
}
names_sigs <- names(gene_sigs_list)
cat("Signatures:", paste(names_sigs, collapse = ", "), "\n")

# --- 3. Source all R functions ---
for (f in list.files("sigQC-master/R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# --- 4. Run the pipeline ---
out_dir <- "results_r"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_all_plots(
  gene_sigs_list  = gene_sigs_list,
  mRNA_expr_matrix = mRNA_expr_matrix,
  names_sigs       = names_sigs,
  names_datasets   = names_datasets,
  out_dir          = out_dir,
  showResults      = FALSE,
  doNegativeControl = FALSE
)

cat("Done. Results in:", out_dir, "\n")
