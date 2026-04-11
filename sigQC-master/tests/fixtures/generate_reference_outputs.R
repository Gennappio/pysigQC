# generate_reference_outputs.R
# Runs all compute_*() functions on the fixture data and saves outputs
# as CSV files for Python cross-validation.
# Run from sigQC-master/: Rscript tests/fixtures/generate_reference_outputs.R

library(mclust)
library(GSVA)
library(biclust)

for (f in list.files("R_refactored", pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

fixture <- readRDS("tests/fixtures/fixture_small.rds")
out_dir <- "tests/fixtures/reference_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- compute_var ----
var_r <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                     fixture$mRNA_expr_matrix, fixture$names_datasets)
for (sig in fixture$names_sigs) {
  for (ds in fixture$names_datasets) {
    rv <- var_r$radar_values[[sig]][[ds]]
    write.csv(t(as.data.frame(rv)),
              file.path(out_dir, paste0("var_radar_", sig, "_", ds, ".csv")),
              row.names = FALSE)
    write.csv(var_r$mean_sd_tables[[sig]][[ds]],
              file.path(out_dir, paste0("var_mean_sd_", sig, "_", ds, ".csv")))
  }
}

# ---- compute_expr ----
expr_r <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                       fixture$mRNA_expr_matrix, fixture$names_datasets)
for (sig in fixture$names_sigs) {
  for (ds in fixture$names_datasets) {
    rv <- expr_r$radar_values[[sig]][[ds]]
    write.csv(t(as.data.frame(rv)),
              file.path(out_dir, paste0("expr_radar_", sig, "_", ds, ".csv")),
              row.names = FALSE)
    write.csv(as.data.frame(expr_r$na_proportions[[sig]][[ds]]),
              file.path(out_dir, paste0("expr_na_prop_", sig, "_", ds, ".csv")))
  }
}
write.csv(data.frame(dataset = names(expr_r$thresholds), threshold = expr_r$thresholds),
          file.path(out_dir, "expr_thresholds.csv"), row.names = FALSE)

# ---- compute_compactness ----
compact_r <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                 fixture$mRNA_expr_matrix, fixture$names_datasets)
for (sig in fixture$names_sigs) {
  for (ds in fixture$names_datasets) {
    rv <- compact_r$radar_values[[sig]][[ds]]
    write.csv(t(as.data.frame(rv)),
              file.path(out_dir, paste0("compact_radar_", sig, "_", ds, ".csv")),
              row.names = FALSE)
    write.csv(compact_r$autocor_matrices[[sig]][[ds]],
              file.path(out_dir, paste0("compact_autocor_", sig, "_", ds, ".csv")))
  }
}

# ---- compute_stan ----
stan_r <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                       fixture$mRNA_expr_matrix, fixture$names_datasets)
for (sig in fixture$names_sigs) {
  for (ds in fixture$names_datasets) {
    rv <- stan_r$radar_values[[sig]][[ds]]
    write.csv(t(as.data.frame(rv)),
              file.path(out_dir, paste0("stan_radar_", sig, "_", ds, ".csv")),
              row.names = FALSE)
    write.csv(data.frame(med_scores = stan_r$med_scores[[sig]][[ds]],
                         z_transf_scores = stan_r$z_transf_scores[[sig]][[ds]]),
              file.path(out_dir, paste0("stan_scores_", sig, "_", ds, ".csv")))
  }
}

# ---- compute_metrics ----
metrics_r <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                             fixture$mRNA_expr_matrix, fixture$names_datasets)
for (sig in fixture$names_sigs) {
  for (ds in fixture$names_datasets) {
    rv <- metrics_r$radar_values[[sig]][[ds]]
    write.csv(t(as.data.frame(rv)),
              file.path(out_dir, paste0("metrics_radar_", sig, "_", ds, ".csv")),
              row.names = FALSE)
    sc <- metrics_r$scores[[sig]][[ds]]
    scores_df <- data.frame(
      med_scores = sc$med_scores,
      mean_scores = sc$mean_scores
    )
    if (!is.null(sc$pca1_scores)) {
      # PCA1 may have fewer rows if NA rows were dropped
      if (length(sc$pca1_scores) == nrow(scores_df)) {
        scores_df$pca1_scores <- sc$pca1_scores
      }
    }
    write.csv(scores_df,
              file.path(out_dir, paste0("metrics_scores_", sig, "_", ds, ".csv")))
  }
}

# ---- compute_radar (full pipeline) ----
radar_values <- list()
for (sig in fixture$names_sigs) {
  radar_values[[sig]] <- list()
  for (ds in fixture$names_datasets) {
    radar_values[[sig]][[ds]] <- c(
      var_r$radar_values[[sig]][[ds]],
      expr_r$radar_values[[sig]][[ds]],
      compact_r$radar_values[[sig]][[ds]],
      metrics_r$radar_values[[sig]][[ds]],
      stan_r$radar_values[[sig]][[ds]]
    )
  }
}
radar_result <- compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)
write.csv(radar_result$output_table,
          file.path(out_dir, "radar_output_table.csv"))
write.csv(data.frame(label = radar_result$legend_labels, area = radar_result$areas),
          file.path(out_dir, "radar_areas.csv"), row.names = FALSE)

cat("Reference outputs generated in", out_dir, "\n")
cat("Files:\n")
cat(paste(" ", list.files(out_dir), collapse = "\n"), "\n")
