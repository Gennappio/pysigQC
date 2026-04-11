# eval_var_loc.R — REFACTORED
#
# Evaluates expression variability of signature genes relative to all genes.
# Split into compute_var() (pure) + plot_var() (side effects) for testability.

# ============================================================================
# compute_var: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $radar_values   — nested list [[sig]][[dataset]] with 6 radar metrics
#   $mean_sd_tables — nested list [[sig]][[dataset]] with mean/SD per gene
#   $all_sd         — nested list [[sig]][[dataset]] with all-gene SD vectors
#   $all_mean       — nested list [[sig]][[dataset]] with all-gene mean vectors
#   $inter          — nested list [[sig]][[dataset]] with intersected gene names
compute_var <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets) {
  radar_values <- list()
  mean_sd_tables <- list()
  all_sd <- list()
  all_mean <- list()
  inter_genes <- list()

  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    radar_values[[names_sigs[k]]] <- list()
    mean_sd_tables[[names_sigs[k]]] <- list()
    all_sd[[names_sigs[k]]] <- list()
    all_mean[[names_sigs[k]]] <- list()
    inter_genes[[names_sigs[k]]] <- list()

    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      inter_genes[[names_sigs[k]]][[names_datasets[i]]] <- inter

      sd_genes <- as.matrix(apply(data.matrix, 1, function(x) stats::sd(as.numeric(x), na.rm = TRUE)))
      mean_genes <- as.matrix(apply(data.matrix, 1, function(x) mean(as.numeric(x), na.rm = TRUE)))

      all_sd[[names_sigs[k]]][[names_datasets[i]]] <- sd_genes
      all_mean[[names_sigs[k]]][[names_datasets[i]]] <- mean_genes

      # Fractional ratio of median SD: sig / (all + sig)
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['sd_median_ratio'] <-
        stats::median(stats::na.omit(sd_genes[inter, 1])) /
        (stats::median(stats::na.omit(sd_genes)) + stats::median(stats::na.omit(sd_genes[inter, 1])))

      # Fractional ratio of absolute skewness
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['abs_skewness_ratio'] <-
        abs(moments::skewness(stats::na.omit(mean_genes[inter, 1]))) /
        (abs(moments::skewness(stats::na.omit(mean_genes))) +
         abs(moments::skewness(stats::na.omit(mean_genes[inter, 1]))))

      # Mean/SD table for this sig-dataset pair
      tbl <- cbind(mean_genes[inter, 1], sd_genes[inter, 1])
      colnames(tbl) <- c("Mean", "SD")
      mean_sd_tables[[names_sigs[k]]][[names_datasets[i]]] <- tbl

      # Coefficient of variation analysis
      coeff_of_var <- as.matrix(apply(data.matrix, 1, function(x) {
        stats::sd(as.numeric(stats::na.omit(x)), na.rm = TRUE) /
          mean(as.numeric(stats::na.omit(x)), na.rm = TRUE)
      }))
      coeff_of_var_gene_sig <- coeff_of_var[inter, 1]
      quantiles_considered <- stats::quantile(coeff_of_var, probs = c(0.9, 0.75, 0.5), na.rm = TRUE)

      radar_values[[names_sigs[k]]][[names_datasets[i]]]['prop_top_10_percent'] <-
        sum(stats::na.omit(coeff_of_var_gene_sig) >= quantiles_considered[1]) /
        length(stats::na.omit(coeff_of_var_gene_sig))
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['prop_top_25_percent'] <-
        sum(stats::na.omit(coeff_of_var_gene_sig) >= quantiles_considered[2]) /
        length(stats::na.omit(coeff_of_var_gene_sig))
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['prop_top_50_percent'] <-
        sum(stats::na.omit(coeff_of_var_gene_sig) >= quantiles_considered[3]) /
        length(stats::na.omit(coeff_of_var_gene_sig))

      radar_values[[names_sigs[k]]][[names_datasets[i]]]['coeff_of_var_ratio'] <-
        abs(abs(stats::median(stats::na.omit(coeff_of_var_gene_sig)))) /
        (abs(stats::median(stats::na.omit(coeff_of_var))) +
         abs(stats::median(stats::na.omit(coeff_of_var_gene_sig))))
    }
  }

  list(
    radar_values = radar_values,
    mean_sd_tables = mean_sd_tables,
    all_sd = all_sd,
    all_mean = all_mean,
    inter = inter_genes
  )
}

# ============================================================================
# plot_var: Side effects — PDFs + CSV files
# ============================================================================
plot_var <- function(compute_result, names_sigs, names_datasets,
                     out_dir, showResults = FALSE) {
  num_rows <- length(names_sigs)
  num_cols <- length(names_datasets)

  max_title_length <- max(nchar(outer(names_datasets, names_sigs, paste)))

  # --- Mean vs SD scatter plots ---
  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_mean_vs_sd.pdf'),
                   width = 4 * num_cols, height = 4 * num_rows)
  }
  graphics::par(mfrow = c(num_rows, num_cols))

  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      sd_genes <- compute_result$all_sd[[names_sigs[k]]][[names_datasets[i]]]
      mean_genes <- compute_result$all_mean[[names_sigs[k]]][[names_datasets[i]]]
      inter <- compute_result$inter[[names_sigs[k]]][[names_datasets[i]]]

      graphics::plot(mean_genes, sd_genes, pch = 19, col = 'grey', cex = 0.5,
                     main = paste0('Mean vs SD of expression\n', names_datasets[i], ' ', names_sigs[k]),
                     xlab = 'Mean', ylab = 'Standard deviation',
                     cex.main = min(1, 4 * 10 / max_title_length))
      graphics::points(mean_genes[inter, 1], sd_genes[inter, 1], pch = 19, col = 'red', cex = 0.5)

      quants_mean <- stats::quantile(mean_genes * is.finite(mean_genes),
                                     probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
      graphics::abline(v = quants_mean[1], lty = 5, col = 'cadetblue1')
      graphics::abline(v = quants_mean[2], lty = 5, col = 'dodgerblue1')
      graphics::abline(v = quants_mean[3], lty = 5, col = 'darkblue')
      graphics::abline(v = quants_mean[4], lty = 5, col = 'dodgerblue1')
      graphics::abline(v = quants_mean[5], lty = 5, col = 'cadetblue1')

      quants_sd <- stats::quantile(sd_genes * is.finite(sd_genes),
                                   probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
      graphics::abline(h = quants_sd[1], lty = 5, col = 'cadetblue1')
      graphics::abline(h = quants_sd[2], lty = 5, col = 'dodgerblue1')
      graphics::abline(h = quants_sd[3], lty = 5, col = 'darkblue')
      graphics::abline(h = quants_sd[4], lty = 5, col = 'dodgerblue1')
      graphics::abline(h = quants_sd[5], lty = 5, col = 'cadetblue1')

      graphics::legend('topright', legend = c('10%, 90%', '25%, 75%', '50%'),
                       col = c('cadetblue1', 'dodgerblue1', 'darkblue'), lty = 5, bty = 'n', cex = 0.5)
    }
  }

  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_mean_vs_sd.pdf'),
                        width = 4 * num_cols, height = 4 * num_rows)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # --- Write mean/SD tables ---
  if (!dir.exists(file.path(out_dir, 'mean_sd_tables'))) {
    dir.create(file.path(out_dir, 'mean_sd_tables'))
  }
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      utils::write.table(
        compute_result$mean_sd_tables[[names_sigs[k]]][[names_datasets[i]]],
        file = file.path(out_dir, 'mean_sd_tables',
                         paste0('mean_sd_table_', names_sigs[k], '_', names_datasets[i], '.txt')),
        quote = FALSE, sep = '\t')
    }
  }
}

# ============================================================================
# eval_var_loc: Original API wrapper — backward-compatible
# ============================================================================
eval_var_loc <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                         out_dir = '~', file = NULL, showResults = FALSE, radar_plot_values) {

  result <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
  plot_var(result, names_sigs, names_datasets, out_dir, showResults)

  cat('Mean vs SD graphs created successfully.\n', file = file)
  cat('Mean vs SD tables written to file successfully.\n', file = file)

  # Merge computed radar values into the accumulator
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
