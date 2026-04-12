# fixture_medium_generator.R
#
# Generates a medium-scale deterministic test fixture for numerical stability
# testing at more realistic dimensions.
#
# Run from sigQC-master/:
#   Rscript tests/fixtures/fixture_medium_generator.R
#
# Dimensions:
#   - 100 genes x 50 samples per dataset
#   - 3 datasets (including one with batch effect)
#   - 5 signatures with different properties:
#       1. "tight_coexpr"   — 8 tightly co-expressed genes
#       2. "loose_coexpr"   — 10 loosely co-expressed genes
#       3. "mixed_quality"  — 6 genes, some expressed, some mostly NA
#       4. "high_var"       — 7 highly variable genes
#       5. "random_baseline"— 12 random genes (negative control)
#
# Edge cases exercised:
#   - Genes shared across signatures (gene overlap)
#   - Genes present in some datasets but not others
#   - Batch effects between datasets
#   - Zero-variance genes, near-zero-variance genes
#   - Varying NA patterns (scattered, block, full-row)
#   - Large signature (12 genes) vs small (6 genes)

set.seed(2024)

n_genes <- 100
n_samples <- 50
gene_names <- paste0("GENE", sprintf("%03d", 1:n_genes))
sample_names_A <- paste0("S_A_", sprintf("%02d", 1:n_samples))
sample_names_B <- paste0("S_B_", sprintf("%02d", 1:n_samples))
sample_names_C <- paste0("S_C_", sprintf("%02d", 1:n_samples))

# ============================================================================
# Dataset A: clean dataset, moderate expression
# ============================================================================
dataset_A <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(dataset_A) <- gene_names
colnames(dataset_A) <- sample_names_A

# Tight co-expression cluster (genes 1-8): high correlation
base_A <- rnorm(n_samples, mean = 7, sd = 2)
for (g in 1:8) {
  dataset_A[g, ] <- base_A + rnorm(n_samples, sd = 0.3)
}

# Loose co-expression cluster (genes 9-18): moderate correlation
base_loose_A <- rnorm(n_samples, mean = 5, sd = 1.5)
for (g in 9:18) {
  dataset_A[g, ] <- base_loose_A + rnorm(n_samples, sd = 1.0)
}

# High-variance genes (genes 19-25)
for (g in 19:25) {
  dataset_A[g, ] <- rnorm(n_samples, mean = runif(1, 3, 10), sd = runif(1, 3, 6))
}

# Remaining genes: random expression
for (g in 26:n_genes) {
  dataset_A[g, ] <- rnorm(n_samples, mean = runif(1, 2, 9), sd = runif(1, 0.5, 2.5))
}

# Inject scattered NAs (5 random positions)
na_positions_A <- cbind(
  sample(1:n_genes, 5),
  sample(1:n_samples, 5)
)
for (i in 1:nrow(na_positions_A)) {
  dataset_A[na_positions_A[i, 1], na_positions_A[i, 2]] <- NA
}

# ============================================================================
# Dataset B: batch-shifted version of A + extra noise
# ============================================================================
dataset_B <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(dataset_B) <- gene_names
colnames(dataset_B) <- sample_names_B

# Same structure but with batch offset (+1.5) and more noise
base_B <- rnorm(n_samples, mean = 8.5, sd = 2)
for (g in 1:8) {
  dataset_B[g, ] <- base_B + rnorm(n_samples, sd = 0.4)
}
base_loose_B <- rnorm(n_samples, mean = 6.5, sd = 1.5)
for (g in 9:18) {
  dataset_B[g, ] <- base_loose_B + rnorm(n_samples, sd = 1.2)
}
for (g in 19:25) {
  dataset_B[g, ] <- rnorm(n_samples, mean = runif(1, 4, 11), sd = runif(1, 3, 6))
}
for (g in 26:n_genes) {
  dataset_B[g, ] <- rnorm(n_samples, mean = runif(1, 3, 10), sd = runif(1, 0.5, 2.5))
}

# Zero-variance gene
dataset_B[20, ] <- rep(5.0, n_samples)

# Near-zero-variance gene (sd ~ 1e-8)
dataset_B[21, ] <- 6.0 + rnorm(n_samples, sd = 1e-8)

# Block of NAs (gene 15, last 10 samples)
dataset_B[15, 41:50] <- NA

# ============================================================================
# Dataset C: reduced gene set (90 genes — genes 91-100 missing)
# ============================================================================
dataset_C <- matrix(0, nrow = 90, ncol = n_samples)
rownames(dataset_C) <- gene_names[1:90]
colnames(dataset_C) <- sample_names_C

base_C <- rnorm(n_samples, mean = 6, sd = 1.8)
for (g in 1:8) {
  dataset_C[g, ] <- base_C + rnorm(n_samples, sd = 0.35)
}
base_loose_C <- rnorm(n_samples, mean = 4.5, sd = 1.5)
for (g in 9:18) {
  dataset_C[g, ] <- base_loose_C + rnorm(n_samples, sd = 0.9)
}
for (g in 19:25) {
  dataset_C[g, ] <- rnorm(n_samples, mean = runif(1, 2, 9), sd = runif(1, 2, 5))
}
for (g in 26:90) {
  dataset_C[g, ] <- rnorm(n_samples, mean = runif(1, 1, 8), sd = runif(1, 0.5, 2))
}
# One fully-NA gene
dataset_C[50, ] <- NA

# ============================================================================
# Signatures
# ============================================================================
# 1. Tight co-expression: genes 1-8
tight_coexpr <- as.matrix(gene_names[1:8])
colnames(tight_coexpr) <- "tight_coexpr"

# 2. Loose co-expression: genes 9-18
loose_coexpr <- as.matrix(gene_names[9:18])
colnames(loose_coexpr) <- "loose_coexpr"

# 3. Mixed quality: 3 expressed + 1 with NAs + 1 zero-var + 1 not in dataset_C
mixed_quality <- as.matrix(c("GENE001", "GENE010", "GENE015", "GENE020", "GENE021", "GENE095"))
colnames(mixed_quality) <- "mixed_quality"

# 4. High variance: genes 19-25
high_var <- as.matrix(gene_names[19:25])
colnames(high_var) <- "high_var"

# 5. Random baseline: 12 genes from the "noise" region
random_baseline <- as.matrix(gene_names[c(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)])
colnames(random_baseline) <- "random_baseline"

# ============================================================================
# Assemble and save
# ============================================================================
gene_sigs_list <- list(
  tight_coexpr = tight_coexpr,
  loose_coexpr = loose_coexpr,
  mixed_quality = mixed_quality,
  high_var = high_var,
  random_baseline = random_baseline
)

mRNA_expr_matrix <- list(
  dataset_A = dataset_A,
  dataset_B = dataset_B,
  dataset_C = dataset_C
)

names_sigs <- names(gene_sigs_list)
names_datasets <- names(mRNA_expr_matrix)

fixture <- list(
  gene_sigs_list = gene_sigs_list,
  mRNA_expr_matrix = mRNA_expr_matrix,
  names_sigs = names_sigs,
  names_datasets = names_datasets
)

script_args <- commandArgs(trailingOnly = FALSE)
file_flag   <- grep("^--file=", script_args, value = TRUE)
if (length(file_flag) > 0) {
  out_dir <- dirname(normalizePath(sub("^--file=", "", file_flag[1])))
} else {
  out_dir <- getwd()
}

# Save RDS
saveRDS(fixture, file = file.path(out_dir, "fixture_medium.rds"))

# Save CSV for Python
for (ds_name in names_datasets) {
  write.csv(mRNA_expr_matrix[[ds_name]],
            file = file.path(out_dir, paste0("fixture_medium_", ds_name, ".csv")))
}

# Save signatures (long format)
sig_rows <- lapply(names(gene_sigs_list), function(sig_name) {
  data.frame(signature = sig_name, gene = gene_sigs_list[[sig_name]][, 1],
             stringsAsFactors = FALSE)
})
sig_df <- do.call(rbind, sig_rows)
write.csv(sig_df, file = file.path(out_dir, "fixture_medium_signatures.csv"),
          row.names = FALSE)

cat("Medium fixture generated:\n")
cat(sprintf("  %d genes, %d samples/dataset, %d datasets, %d signatures\n",
            n_genes, n_samples, length(names_datasets), length(names_sigs)))
cat(sprintf("  Edge cases: NAs (scattered, block, full-row), zero-var, near-zero-var,\n"))
cat(sprintf("              missing genes in dataset_C, gene overlap via mixed_quality\n"))
