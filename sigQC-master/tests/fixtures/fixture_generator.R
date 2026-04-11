# fixture_generator.R
#
# Generates deterministic test fixtures for sigQC unit tests.
# Run this script to regenerate fixture_small.rds and fixture_small.csv.
#
# The fixtures are designed to exercise edge cases:
#   - Signature 1 ("compact_sig"): 5 genes with high inter-gene correlation
#   - Signature 2 ("diffuse_sig"): 5 genes with low/no correlation
#   - Dataset 1 ("dataset_A"): 10 samples, includes 1 NA value
#   - Dataset 2 ("dataset_B"): 10 samples, includes 1 zero-variance gene
#   - Gene "gene_missing" appears in compact_sig but not in any dataset (tests padding)

set.seed(42)

n_samples <- 10
sample_names <- paste0("sample_", seq_len(n_samples))

# --- Dataset A: 20 genes x 10 samples ---
# First 4 genes are correlated (for compact signature), rest are random
base_signal <- rnorm(n_samples, mean = 5, sd = 2)
dataset_A <- matrix(0, nrow = 20, ncol = n_samples)
rownames(dataset_A) <- paste0("gene_", 1:20)
colnames(dataset_A) <- sample_names

# Compact signature genes (1-4): correlated with base_signal + noise
for (g in 1:4) {
  dataset_A[g, ] <- base_signal + rnorm(n_samples, sd = 0.5)
}
# Random genes (5-20)
for (g in 5:20) {
  dataset_A[g, ] <- rnorm(n_samples, mean = runif(1, 2, 8), sd = runif(1, 0.5, 3))
}
# Inject one NA value
dataset_A[3, 5] <- NA

# --- Dataset B: 20 genes x 10 samples ---
dataset_B <- matrix(0, nrow = 20, ncol = n_samples)
rownames(dataset_B) <- paste0("gene_", 1:20)
colnames(dataset_B) <- sample_names

base_signal_B <- rnorm(n_samples, mean = 6, sd = 1.5)
for (g in 1:4) {
  dataset_B[g, ] <- base_signal_B + rnorm(n_samples, sd = 0.3)
}
for (g in 5:20) {
  dataset_B[g, ] <- rnorm(n_samples, mean = runif(1, 3, 7), sd = runif(1, 0.3, 2))
}
# Make gene_5 zero-variance (constant expression) to test BUG-2 fix
dataset_B[5, ] <- rep(4.0, n_samples)

# --- Gene signatures ---
# Compact: 4 correlated genes + 1 missing gene
compact_sig <- as.matrix(c("gene_1", "gene_2", "gene_3", "gene_4", "gene_missing"))
colnames(compact_sig) <- "compact_sig"

# Diffuse: 5 random/uncorrelated genes
diffuse_sig <- as.matrix(c("gene_5", "gene_10", "gene_12", "gene_15", "gene_18"))
colnames(diffuse_sig) <- "diffuse_sig"

# --- Assemble fixture ---
gene_sigs_list <- list(
  compact_sig = compact_sig,
  diffuse_sig = diffuse_sig
)

mRNA_expr_matrix <- list(
  dataset_A = dataset_A,
  dataset_B = dataset_B
)

names_sigs <- names(gene_sigs_list)
names_datasets <- names(mRNA_expr_matrix)

fixture <- list(
  gene_sigs_list = gene_sigs_list,
  mRNA_expr_matrix = mRNA_expr_matrix,
  names_sigs = names_sigs,
  names_datasets = names_datasets
)

# Save as RDS
saveRDS(fixture, file = file.path(dirname(sys.frame(1)$ofile %||% "."), "fixture_small.rds"))

# Save expression matrices as CSV (for Python consumption)
for (ds_name in names_datasets) {
  write.csv(mRNA_expr_matrix[[ds_name]],
            file = file.path(dirname(sys.frame(1)$ofile %||% "."),
                             paste0("fixture_", ds_name, ".csv")))
}

# Save signatures as CSV
sig_df <- data.frame(
  signature = c(rep("compact_sig", 5), rep("diffuse_sig", 5)),
  gene = c(compact_sig[, 1], diffuse_sig[, 1])
)
write.csv(sig_df, file = file.path(dirname(sys.frame(1)$ofile %||% "."),
                                    "fixture_signatures.csv"),
          row.names = FALSE)

cat("Fixtures generated successfully.\n")
