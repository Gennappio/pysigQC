# eval_compactness_loc.R — REFACTORED
#
# Evaluates compactness (internal coherence) of gene signatures via Spearman
# autocorrelation matrices. Split into compute_compactness() + plot_compactness().

# ============================================================================
# compute_compactness: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $radar_values       — nested list [[sig]][[dataset]] with autocor_median
#   $autocor_matrices   — nested list [[sig]][[dataset]] with correlation matrices
#   $rank_product_tables — nested list [[sig]] with rank product results (if >1 dataset)
compute_compactness <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix,
                                names_datasets, logged = TRUE, origin = NULL) {
  radar_values <- list()
  autocor_matrices <- list()
  rank_product_tables <- list()

  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    radar_values[[names_sigs[k]]] <- list()
    autocor_matrices[[names_sigs[k]]] <- list()

    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      autocors <- stats::cor(t(stats::na.omit(data.matrix[inter, ])), method = 'spearman')
      autocor_matrices[[names_sigs[k]]][[names_datasets[i]]] <- autocors

      if (dim(autocors)[1] > 1) {
        radar_values[[names_sigs[k]]][[names_datasets[i]]]['autocor_median'] <-
          stats::median(autocors, na.rm = TRUE)
      } else {
        radar_values[[names_sigs[k]]][[names_datasets[i]]]['autocor_median'] <- 0
      }
    }
  }

  # Rank product analysis (only if >1 dataset and RankProd installed)
  if (length(names_datasets) > 1) {
    RankProdInstalled <- (nchar(system.file(package = 'RankProd')) > 0)
    if (RankProdInstalled) {
      for (k in seq_along(names_sigs)) {
        gene_sig <- gene_sigs_list[[names_sigs[k]]]
        if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

        overall_rank_mat <- matrix(NA, nrow = length(unique(gene_sig)), ncol = length(names_datasets))
        row.names(overall_rank_mat) <- unique(gene_sig)
        colnames(overall_rank_mat) <- names_datasets

        for (i in seq_along(names_datasets)) {
          data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
          inter <- intersect(unique(gene_sig), rownames(data.matrix))
          autocors <- stats::cor(t(stats::na.omit(data.matrix[inter, ])), method = 'spearman')
          median_scores <- as.matrix(apply(autocors, 2, function(x) stats::median(stats::na.omit(x))))
          overall_rank_mat[rownames(median_scores), i] <- median_scores[, 1]
        }

        if (is.null(origin)) origin <- rep(1, length(names_datasets))
        RP.out <- RankProd::RPadvance(data = overall_rank_mat,
                                       cl = rep(1, times = length(names_datasets)),
                                       origin = origin, logged = logged,
                                       gene.names = rownames(overall_rank_mat))
        table_rp <- cbind(RP.out$pfp, RP.out$pval, RP.out$RPs)
        colnames(table_rp) <- c(paste0("pfp_", colnames(RP.out$pfp)),
                                 paste0('p_val', colnames(RP.out$pval)),
                                 paste0("Rank_Product_", colnames(RP.out$RPs)))
        rank_product_tables[[names_sigs[k]]] <- list(RP.out = RP.out, table = table_rp)
      }
    }
  }

  list(
    radar_values = radar_values,
    autocor_matrices = autocor_matrices,
    rank_product_tables = rank_product_tables
  )
}

# ============================================================================
# plot_compactness: Side effects — PDFs + CSV files
# ============================================================================
plot_compactness <- function(compute_result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
                             names_datasets, out_dir, file = NULL, showResults = FALSE) {
  max_title_length <- max(nchar(outer(names_datasets, names_sigs, paste)))

  # --- Autocorrelation heatmaps ---
  grDevices::dev.new()
  graphics::par(cex.main = min(0.8, (3 * 6 / max_title_length)), cex.lab = 0.6,
                oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0))

  hmaps <- lapply(seq_len(length(names_sigs) * length(names_datasets)), function(idx) {
    dataset_ind <- ((idx - 1) %% length(names_datasets)) + 1
    sig_ind <- ceiling(idx / length(names_datasets))
    autocors <- compute_result$autocor_matrices[[names_sigs[sig_ind]]][[names_datasets[dataset_ind]]]

    # Write autocorrelation matrix
    if (!dir.exists(file.path(out_dir, 'autocorrelation_matrices'))) {
      dir.create(file.path(out_dir, 'autocorrelation_matrices'))
    }
    utils::write.table(autocors,
                       file = file.path(out_dir, 'autocorrelation_matrices',
                                        paste0('autocorrelation_matrix_', names_sigs[sig_ind], '_',
                                               names_datasets[dataset_ind], '.txt')),
                       quote = FALSE, sep = '\t')

    tryCatch({
      gplots::heatmap.2(stats::na.omit(autocors),
                         col = gplots::colorpanel(100, "blue", "white", "red"),
                         trace = "none", na.color = "grey",
                         labRow = rownames(autocors), labCol = colnames(autocors),
                         main = paste0("\n\nIntra-sig. Corr.\n",
                                       names_datasets[[dataset_ind]], ' ', names_sigs[[sig_ind]]),
                         dendrogram = "col", symbreaks = TRUE, Rowv = TRUE, Colv = TRUE,
                         key.xlab = 'Rho', key.ylab = NA, key.title = NA,
                         cexRow = max(min(0.5, (4 * 4 / length(rownames(autocors)))), 0.06),
                         cexCol = max(min(0.5, (4 * 4 / length(rownames(autocors)))), 0.06),
                         margins = c(1 + (max(nchar(rownames(autocors))) / 2),
                                    1 + (max(nchar(rownames(autocors))) / 2)))
    }, error = function(err) {
      graphics::plot.new()
      graphics::title(paste0('\n\nToo many NA values in \n',
                             names_datasets[dataset_ind], ' ', names_sigs[sig_ind]))
      cat(paste0("There was an error,  ", names_datasets[dataset_ind], " ",
                 names_sigs[sig_ind], " ", err, '\n'), file = file)
    })
    grab_grob()
  })

  draw.heatmaps(hmaps, names_datasets, names_sigs)
  grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_autocor_hmps.pdf'),
                      width = 4 * length(names_datasets), height = 4 * length(names_sigs))
  if (grDevices::dev.cur() != 1) grDevices::dev.off()
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # --- Autocorrelation density plots ---
  legend_names <- c()
  legend_cols <- c()
  legend_lty <- c()
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      legend_names <- c(legend_names, paste0(names_datasets[i], ' ', names_sigs[k]))
      legend_cols <- c(legend_cols, i)
      legend_lty <- c(legend_lty, k)
    }
  }

  graphics::par(mar = c(0, 0, 0, 0), cex = 0.6)
  graphics::plot.new()
  l <- graphics::legend(0, 0, legend_names, col = legend_cols, lty = legend_lty,
                        lwd = rep(1, times = (length(names_datasets) * length(names_sigs))),
                        pt.cex = 1, cex = min(0.5, (4 * 10 / max_title_length)), plot = FALSE)
  w <- graphics::grconvertX(l$rect$w, to = 'ndc') - graphics::grconvertX(0, to = 'ndc')
  w <- graphics::grconvertX(w, from = "ndc", to = "inches") +
    graphics::grconvertX(10, from = "device", to = "inches")

  if (showResults) {
    grDevices::dev.new()
  } else {
    grDevices::pdf(file.path(out_dir, 'sig_autocor_dens.pdf'), width = 5, height = 5)
  }
  graphics::par(cex.main = 0.8, cex.lab = 0.6, mar = c(3, 3, 4, 1),
                mfrow = c(1, 1), xpd = TRUE, omi = c(0, 0, 0, w))

  max_dens <- -9999
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      autocors <- compute_result$autocor_matrices[[names_sigs[k]]][[names_datasets[i]]]
      if (dim(autocors)[1] > 1) {
        cur_max <- max(stats::density(unlist(stats::na.omit(autocors)))$y)
        if (max_dens < cur_max) max_dens <- cur_max
      }
    }
  }

  plots_count <- 0
  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      autocors <- compute_result$autocor_matrices[[names_sigs[k]]][[names_datasets[i]]]
      if (dim(autocors)[1] > 1) {
        if (plots_count == 0) {
          graphics::plot(stats::density(unlist(stats::na.omit(autocors))),
                         ylim = c(0, ceiling(max_dens)), xlim = c(-1, 1),
                         col = i, main = NA, lwd = 2, lty = k)
          plots_count <- 1
        } else {
          graphics::lines(stats::density(unlist(stats::na.omit(autocors))),
                          ylim = c(0, ceiling(max_dens)), xlim = c(-1, 1),
                          col = i, main = NA, lwd = 2, lty = k)
        }
      }
    }
  }

  graphics::mtext(side = 2, line = 2, 'Density', cex = 0.8)
  graphics::mtext(side = 1, line = 2, 'Rho', cex = 0.8)
  graphics::mtext(side = 3, line = 2, 'Intra-sig. Corr. Density')
  graphics::par(cex = 0.6)
  graphics::legend(graphics::par('usr')[2] + 0.05, graphics::par('usr')[4], xpd = NA,
                   legend_names, col = legend_cols, lty = legend_lty,
                   lwd = rep(1, times = (length(names_datasets) * length(names_sigs))),
                   pt.cex = 1, cex = min(0.5, (4 * 10 / max_title_length)))

  if (showResults) {
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_autocor_dens.pdf'),
                        width = 5, height = 5)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()

  # --- Rank product plots ---
  if (length(compute_result$rank_product_tables) > 0) {
    for (k in seq_along(names_sigs)) {
      if (!is.null(compute_result$rank_product_tables[[names_sigs[k]]])) {
        if (showResults) {
          grDevices::dev.new()
        } else {
          grDevices::pdf(file.path(out_dir, paste0('sig_autocor_rankProd_', names_sigs[k], '.pdf')),
                         width = 10, height = 10)
        }
        graphics::par(cex.main = 0.8, cex.lab = 0.6, oma = c(2, 2, 2, 2), mar = c(4, 4, 4, 4))
        RankProd::plotRP(compute_result$rank_product_tables[[names_sigs[k]]]$RP.out, cutoff = 0.05)

        if (!dir.exists(file.path(out_dir, 'rank_prod'))) {
          dir.create(file.path(out_dir, 'rank_prod'))
        }
        tbl <- compute_result$rank_product_tables[[names_sigs[k]]]$table
        utils::write.csv(tbl[order(tbl[, 1]), ],
                         file = file.path(out_dir, 'rank_prod',
                                          paste0('rank_product_table1_', names_sigs[k], '.txt')),
                         quote = FALSE)
        utils::write.csv(tbl[order(tbl[, 2]), ],
                         file = file.path(out_dir, 'rank_prod',
                                          paste0('rank_product_table2_', names_sigs[k], '.txt')),
                         quote = FALSE)

        if (showResults) {
          grDevices::dev.copy(grDevices::pdf,
                              file.path(out_dir, paste0('sig_autocor_rankProd_', names_sigs[k], '.pdf')),
                              width = 10, height = 10)
        }
        if (grDevices::dev.cur() != 1) grDevices::dev.off()
      }
    }
  }
}

# ============================================================================
# eval_compactness_loc: Original API wrapper — backward-compatible
# ============================================================================
eval_compactness_loc <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                                  out_dir = '~', file = NULL, showResults = FALSE,
                                  radar_plot_values, logged = TRUE, origin = NULL) {

  result <- compute_compactness(gene_sigs_list, names_sigs, mRNA_expr_matrix,
                                 names_datasets, logged, origin)
  plot_compactness(result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
                   names_datasets, out_dir, file, showResults)

  cat("Intra-sig. corr. metrics successfully computed.\n", file = file)

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
