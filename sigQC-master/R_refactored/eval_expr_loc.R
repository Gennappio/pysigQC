# eval_expr_loc.R — REFACTORED
#
# Evaluates expression level properties of gene signatures across datasets.
# Split into compute_expr() (pure) + plot_expr() (side effects).

# ============================================================================
# compute_expr: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $radar_values     — nested list [[sig]][[dataset]] with 2 radar metrics
#   $na_proportions   — nested list [[sig]][[dataset]] with per-gene NA proportions
#   $expr_proportions — nested list [[sig]][[dataset]] with per-gene expression proportions
#   $thresholds       — named vector of expression thresholds per dataset
compute_expr <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                         thresholds = NULL) {
  radar_values <- list()
  na_proportions <- list()
  expr_proportions <- list()

  # --- NA proportion analysis ---
  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    radar_values[[names_sigs[k]]] <- list()
    na_proportions[[names_sigs[k]]] <- list()
    expr_proportions[[names_sigs[k]]] <- list()

    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      genes_expr <- data.matrix[inter, ]
      gene_expr_vals <- (rowSums(is.na(genes_expr)) / dim(genes_expr)[2])
      gene_expr_vals[setdiff(gene_sig, row.names(data.matrix))] <- 1
      gene_expr_vals <- -sort(-gene_expr_vals)
      na_proportions[[names_sigs[k]]][[names_datasets[i]]] <- gene_expr_vals
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['med_prop_na'] <- stats::median(1 - gene_expr_vals)
    }
  }

  # --- Compute thresholds ---
  if (length(thresholds) == 0) {
    thresholds <- rep(0, length(names_datasets))
    for (i in seq_along(names_datasets)) {
      thresholds[i] <- stats::median(unlist(stats::na.omit(mRNA_expr_matrix[[names_datasets[i]]])))
    }
  }
  if (length(names(thresholds)) == 0) {
    names(thresholds) <- names_datasets
  }

  # --- Expression proportion analysis ---
  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      genes_expr <- data.matrix[inter, ]
      gene_expr_vals <- 1 - ((rowSums(genes_expr < thresholds[i])) / (dim(genes_expr)[2]))
      gene_expr_vals[setdiff(gene_sig, row.names(data.matrix))] <- 0
      gene_expr_vals <- sort(gene_expr_vals)
      expr_proportions[[names_sigs[k]]][[names_datasets[i]]] <- gene_expr_vals
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['med_prop_above_med'] <- stats::median(gene_expr_vals)
    }
  }

  list(
    radar_values = radar_values,
    na_proportions = na_proportions,
    expr_proportions = expr_proportions,
    thresholds = thresholds
  )
}

# ============================================================================
# plot_expr: Side effects — PDFs + CSV files
# ============================================================================
plot_expr <- function(compute_result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
                      names_datasets, out_dir, showResults = FALSE) {
  num_rows <- length(names_sigs)
  num_cols <- length(names_datasets)
  max_line_length <- max(nchar(outer(names_datasets, names_sigs, paste)))
  thresholds <- compute_result$thresholds

  # --- Plot 1: NA proportion barcharts ---
  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_expr_barcharts_NA_values.pdf'),
                   width = 4 * num_cols, height = 4 * num_rows)
  }
  graphics::par(mfrow = c(num_rows, num_cols), cex = 0.7, cex.axis = 0.5)

  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      gene_expr_vals <- compute_result$na_proportions[[names_sigs[k]]][[names_datasets[i]]]
      if (max(gene_expr_vals) == 0) {
        bar_expr <- graphics::barplot(gene_expr_vals, xlab = "Signature Gene IDs",
                                      ylab = "Proportion of NA expression",
                                      main = paste0("Signature gene expression\n", names_datasets[i], ' ', names_sigs[k]),
                                      axisnames = FALSE, axis = FALSE, ylim = c(0, 1),
                                      cex.main = min(1, 4 * 12 / max_line_length), border = NA)
      } else {
        bar_expr <- graphics::barplot(gene_expr_vals, xlab = "Signature Gene IDs",
                                      ylab = "Proportion of NA expression",
                                      main = paste0("Signature gene expression\n", names_datasets[i], ' ', names_sigs[k]),
                                      axisnames = FALSE, axis = FALSE,
                                      cex.main = min(1, 4 * 12 / max_line_length), border = NA)
      }
      graphics::text(bar_expr, graphics::par("usr")[3], labels = names(gene_expr_vals),
                     srt = 45, adj = c(1.1, 1.1), xpd = TRUE,
                     cex = max(min(0.5, (0.5 * 4 * 12) / (sqrt(2) * length(gene_expr_vals))), 0.06))
      graphics::axis(2)
    }
  }
  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_expr_barcharts_NA_values.pdf'),
                        width = 4 * num_cols, height = 4 * num_rows)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # --- Plot 2: Expression proportion barcharts ---
  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_expr_barcharts.pdf'),
                   width = 4 * num_cols, height = 4 * num_rows)
  }
  graphics::par(mfrow = c(num_rows, num_cols), cex = 0.7, cex.axis = 0.5)

  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      gene_expr_vals <- compute_result$expr_proportions[[names_sigs[k]]][[names_datasets[i]]]
      bar_expr <- graphics::barplot(gene_expr_vals, xlab = "Signature Gene IDs",
                                    ylab = "Proportion with expression above threshold",
                                    main = paste0("Signature gene expression\n", names_datasets[i], ' ', names_sigs[k]),
                                    axisnames = FALSE, axis = FALSE,
                                    cex.main = min(1, 4 * 12 / max_line_length), border = NA)
      graphics::text(bar_expr, graphics::par("usr")[3], labels = names(gene_expr_vals),
                     srt = 45, adj = c(1.1, 1.1), xpd = TRUE,
                     cex = max(min(0.5, (0.5 * 4 * 12) / (sqrt(2) * length(gene_expr_vals))), 0.06))
      graphics::axis(2)
    }
  }
  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_expr_barcharts.pdf'),
                        width = 4 * num_cols, height = 4 * num_rows)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # --- Plot 3: Expression density plots ---
  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_expr_density_plots.pdf'),
                   width = 4 * num_cols, height = 4 * num_rows)
  }
  graphics::par(mfrow = c(num_rows, num_cols))

  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      genes_expr <- data.matrix[inter, ]
      gene_expr_vals <- 1 - (rowSums(genes_expr < thresholds[i]) / (dim(genes_expr)[2]))
      gene_expr_vals[setdiff(gene_sig, row.names(data.matrix))] <- 0
      graphics::plot(stats::density(stats::na.omit(gene_expr_vals), adjust = 0.25),
                     main = paste0("Signature gene expression\n", names_datasets[i], ' ', names_sigs[k]),
                     ylab = "Density", cex.main = min(1, 3.5 * 10 / max_line_length))
    }
  }
  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_expr_density_plots.pdf'),
                        width = 4 * num_cols, height = 4 * num_rows)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # --- Write expression tables ---
  if (!dir.exists(file.path(out_dir, 'expression_tables'))) {
    dir.create(file.path(out_dir, 'expression_tables'))
  }
  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      genes_expr <- data.matrix[inter, ]
      gene_expr_vals <- 1 - (rowSums(genes_expr < thresholds[i]) / (dim(genes_expr)[2]))
      gene_expr_vals[setdiff(gene_sig, row.names(data.matrix))] <- 0
      gene_expr_vals_na <- (rowSums(is.na(genes_expr)) / dim(genes_expr)[2])
      gene_expr_vals_na[setdiff(gene_sig, row.names(data.matrix))] <- 1
      output_table <- cbind(gene_expr_vals, gene_expr_vals_na)
      colnames(output_table) <- c('Proportion_above_threshold', 'NA proportion')
      utils::write.table(output_table,
                         file = file.path(out_dir, 'expression_tables',
                                          paste0('expression_table_', names_sigs[k], '_', names_datasets[i], '.txt')),
                         quote = FALSE, sep = '\t')
    }
  }
}

# ============================================================================
# eval_expr_loc: Original API wrapper — backward-compatible
# ============================================================================
eval_expr_loc <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                          thresholds = NULL, out_dir = '~', file = NULL,
                          showResults = FALSE, radar_plot_values) {

  result <- compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, thresholds)
  plot_expr(result, gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, out_dir, showResults)

  cat('Expression and density graphs created successfully.\n', file = file)

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
