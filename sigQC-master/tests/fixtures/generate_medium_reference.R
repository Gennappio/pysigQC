# generate_medium_reference.R
# Generates reference outputs from the medium fixture for Python cross-validation.
# Run from sigQC-master/: Rscript tests/fixtures/generate_medium_reference.R

library(mclust)
library(GSVA)
library(biclust)

for (f in list.files("R_refactored", pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

fixture <- readRDS("tests/fixtures/fixture_medium.rds")
out_dir <- "tests/fixtures/reference_outputs_medium"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ns <- fixture$names_sigs
nd <- fixture$names_datasets

# Run all compute modules
cat("Running compute_var...\n")
var_r <- compute_var(fixture$gene_sigs_list, ns, fixture$mRNA_expr_matrix, nd)

cat("Running compute_expr...\n")
expr_r <- compute_expr(fixture$gene_sigs_list, ns, fixture$mRNA_expr_matrix, nd)

cat("Running compute_compactness...\n")
compact_r <- compute_compactness(fixture$gene_sigs_list, ns, fixture$mRNA_expr_matrix, nd)

cat("Running compute_stan...\n")
stan_r <- compute_stan(fixture$gene_sigs_list, ns, fixture$mRNA_expr_matrix, nd)

cat("Running compute_metrics...\n")
metrics_r <- compute_metrics(fixture$gene_sigs_list, ns, fixture$mRNA_expr_matrix, nd)

# Save radar values per module
for (sig in ns) {
  for (ds in nd) {
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

# Assemble full radar and save output table
cat("Assembling radar chart...\n")
radar_values <- list()
for (sig in ns) {
  radar_values[[sig]] <- list()
  for (ds in nd) {
    radar_values[[sig]][[ds]] <- c(
      var_r$radar_values[[sig]][[ds]],
      expr_r$radar_values[[sig]][[ds]],
      compact_r$radar_values[[sig]][[ds]],
      metrics_r$radar_values[[sig]][[ds]],
      stan_r$radar_values[[sig]][[ds]]
    )
  }
}
radar_result <- compute_radar(radar_values, ns, nd)
write.csv(radar_result$output_table,
          file.path(out_dir, "radar_output_table.csv"))

cat("Medium reference outputs saved to", out_dir, "\n")
cat(length(list.files(out_dir)), "files generated.\n")
