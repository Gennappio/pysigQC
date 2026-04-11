# eval_stan_loc.R — REFACTORED
#
# Evaluates the effect of z-score standardization on signature scoring.
# Split into compute_stan() (pure) + plot_stan() (side effects).

# ============================================================================
# compute_stan: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $radar_values     — nested list [[sig]][[dataset]] with standardization_comp
#   $med_scores       — nested list [[sig]][[dataset]] with raw median scores
#   $z_transf_scores  — nested list [[sig]][[dataset]] with z-median scores
compute_stan <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets) {
  radar_values <- list()
  med_scores_all <- list()
  z_transf_scores_all <- list()

  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    radar_values[[names_sigs[k]]] <- list()
    med_scores_all[[names_sigs[k]]] <- list()
    z_transf_scores_all[[names_sigs[k]]] <- list()

    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))

      z_transf_mRNA <- data.matrix[inter, ]
      for (gene in inter) {
        gene_values <- as.numeric(z_transf_mRNA[gene, ])
        gene_sd <- stats::sd(gene_values, na.rm = TRUE)
        if (is.na(gene_sd) || gene_sd == 0) {
          z_transf_mRNA[gene, ] <- 0
        } else {
          z_transf_mRNA[gene, ] <- (gene_values - mean(gene_values, na.rm = TRUE)) / gene_sd
        }
      }

      z_transf_scores <- apply(z_transf_mRNA[inter, ], 2, function(x) stats::median(stats::na.omit(x)))
      med_scores <- apply(data.matrix[inter, ], 2, function(x) stats::median(stats::na.omit(x)))

      rho <- stats::cor(med_scores, z_transf_scores, method = 'spearman')

      radar_values[[names_sigs[k]]][[names_datasets[i]]]['standardization_comp'] <- rho
      med_scores_all[[names_sigs[k]]][[names_datasets[i]]] <- med_scores
      z_transf_scores_all[[names_sigs[k]]][[names_datasets[i]]] <- z_transf_scores
    }
  }

  list(
    radar_values = radar_values,
    med_scores = med_scores_all,
    z_transf_scores = z_transf_scores_all
  )
}

# ============================================================================
# plot_stan: Side effects — PDFs + CSV files
# ============================================================================
plot_stan <- function(compute_result, names_sigs, names_datasets,
                      out_dir, showResults = FALSE) {
  num_rows <- length(names_sigs)
  num_cols <- length(names_datasets)
  max_title_length <- max(nchar(outer(names_datasets, names_sigs, paste)))

  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_standardisation_comp.pdf'),
                   width = 4 * num_cols, height = 4 * num_rows)
  }
  graphics::par(mfrow = c(num_rows, num_cols), oma = c(2, 2, 2, 2), mar = c(4, 4, 4, 4))

  if (!dir.exists(file.path(out_dir, 'standardisation_tables'))) {
    dir.create(file.path(out_dir, 'standardisation_tables'))
  }

  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      med_scores <- compute_result$med_scores[[names_sigs[k]]][[names_datasets[i]]]
      z_transf_scores <- compute_result$z_transf_scores[[names_sigs[k]]][[names_datasets[i]]]
      rho <- compute_result$radar_values[[names_sigs[k]]][[names_datasets[i]]]['standardization_comp']

      jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      graphics::smoothScatter(med_scores, z_transf_scores, colramp = jet.colors,
                              xlab = NA, ylab = NA, main = 'Median vs Z-median')
      graphics::points(med_scores, z_transf_scores, pch = '.')
      graphics::par(new = TRUE, srt = 45)
      graphics::plot(stats::density(med_scores), axes = FALSE, xlab = NA, ylab = NA,
                     col = 'red', main = NA, lwd = 2)
      graphics::axis(side = 4)
      graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
      graphics::mtext(side = 2, line = 2, 'Z-median', cex = 0.8)
      graphics::mtext(side = 1, line = 2, 'Median', cex = 0.8)
      graphics::mtext(side = 3, line = 2.5,
                      paste0(names_datasets[i], ', ', names_sigs[k]),
                      cex = min(1, 4 * 10 / max_title_length))
      graphics::mtext(paste0('rho = ', format(rho, digits = 2)),
                      side = 3, line = 0, cex = 0.6, at = max(med_scores))

      output_mat <- cbind(med_scores, z_transf_scores)
      colnames(output_mat) <- c('Median Scores', 'Z-Median Scores')
      utils::write.table(output_mat,
                         file = file.path(out_dir, 'standardisation_tables',
                                          paste0('standardisation_table_', names_sigs[k], '_', names_datasets[i], '.txt')),
                         quote = FALSE, sep = '\t')
    }
  }

  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_standardisation_comp.pdf'),
                        width = 4 * num_cols, height = 4 * num_rows)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()
}

# ============================================================================
# eval_stan_loc: Original API wrapper — backward-compatible
# ============================================================================
eval_stan_loc <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                          out_dir = '~', file = NULL, showResults = FALSE, radar_plot_values) {

  result <- compute_stan(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
  plot_stan(result, names_sigs, names_datasets, out_dir, showResults)

  cat('Standardisation compared successfully.\n', file = file)

  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      for (metric in names(result$radar_values[[names_sigs[k]]][[names_datasets[i]]])) {
        radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][metric] <-
          result$radar_values[[names_sigs[k]]][[names_datasets[i]]][metric]
      }
    }
  }
  radar_plot_values
}
