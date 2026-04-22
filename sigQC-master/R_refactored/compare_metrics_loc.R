# compare_metrics_loc.R — REFACTORED
#
# Compares various signature summary scoring metrics (Mean, Median, PCA1, GSVA,
# ssGSEA, PLAGE). Split into compute_metrics() + plot_metrics().
# 5 sections: (1) Mean/Median/PCA1 correlations, (2) ES comparisons,
# (3) Correlation heatmaps, (4) Gaussian mixture models, (5) QQ plots.

# ============================================================================
# compute_metrics: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $radar_values    — nested list [[sig]][[dataset]] with 4 radar metrics
#   $scores          — nested list [[sig]][[dataset]] with all score vectors
#   $pca_results     — nested list [[sig]][[dataset]] with PCA objects
#   $score_cor_mats  — named list of scoring metric correlation matrices
#   $mixture_models  — nested list [[sig]][[dataset]] with mclust results
#
# When `radar_only = TRUE`, the caller only needs the 4 radar metrics; the
# GSVA/ssGSEA/PLAGE enrichment scores, the score correlation matrices, and
# the Gaussian mixture models are skipped. This is used by the negative /
# permutation control pipelines where these plotting-only outputs are never
# consumed, and reduces the cost of each resampling iteration by >20x.
compute_metrics <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix,
                            names_datasets, radar_only = FALSE) {
  radar_values <- list()
  scores_all <- list()
  pca_results <- list()
  score_cor_mats <- list()
  mixture_models <- list()

  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    radar_values[[names_sigs[k]]] <- list()
    scores_all[[names_sigs[k]]] <- list()
    pca_results[[names_sigs[k]]] <- list()
    mixture_models[[names_sigs[k]]] <- list()

    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      data.matrix[!(is.finite(as.matrix(data.matrix)))] <- NA
      inter <- intersect(gene_sig, rownames(data.matrix))

      # --- Compute summary scores ---
      med_scores <- apply(data.matrix[inter, ], 2, function(x) stats::median(stats::na.omit(x)))
      mean_scores <- apply(data.matrix[inter, ], 2, function(x) mean(stats::na.omit(x)))

      pca1_scores <- NULL
      pca_obj <- NULL
      props_of_variances <- NULL
      tryCatch({
        pca_obj <- stats::prcomp(t(stats::na.omit(data.matrix[inter, ])), retx = TRUE)
        props_of_variances <- pca_obj$sdev^2 / (sum(pca_obj$sdev^2))
        pca1_scores <- pca_obj$x[, 1]
      }, error = function(e) {
        pca1_scores <<- NULL
        pca_obj <<- NULL
      })

      pca_results[[names_sigs[k]]][[names_datasets[i]]] <- list(
        pca_obj = pca_obj,
        props_of_variances = props_of_variances
      )

      # --- Compute enrichment scores ---
      es.gsva <- es.ssGSEA <- es.plage <- NULL
      if (!radar_only) {
        data.matrix.gsva <- data.matrix
        if (!is.matrix(data.matrix)) {
          data.matrix.gsva <- as.matrix.data.frame(data.matrix, rownames.force = TRUE)
        }
        tryCatch({
          gsvaPar <- GSVA::ssgseaParam(data.matrix.gsva, list(inter))
          es.ssGSEA <- suppressWarnings(GSVA::gsva(gsvaPar, verbose = FALSE))[1, ]

          gsvaPar <- GSVA::gsvaParam(data.matrix.gsva, list(inter))
          es.gsva <- suppressWarnings(GSVA::gsva(gsvaPar, verbose = FALSE))[1, ]

          gsvaPar <- GSVA::plageParam(data.matrix.gsva, list(inter))
          es.plage <- suppressWarnings(GSVA::gsva(gsvaPar, verbose = FALSE))[1, ]
        }, error = function(e) {
          # GSVA may fail with too few genes
        })
      }

      # --- Compute correlations ---
      if (!is.null(pca1_scores)) {
        common_score_cols <- intersect(names(med_scores),
                                        intersect(names(mean_scores), names(pca1_scores)))
      } else {
        common_score_cols <- intersect(names(med_scores), names(mean_scores))
      }

      if (length(common_score_cols) > 1) {
        rho_mean_med <- stats::cor(med_scores[common_score_cols],
                                    mean_scores[common_score_cols], method = 'spearman')
      } else {
        rho_mean_med <- 0
      }

      if (length(pca1_scores) > 1) {
        rho_mean_pca1 <- stats::cor(mean_scores[common_score_cols],
                                     pca1_scores[common_score_cols], method = 'spearman')
        rho_pca1_med <- stats::cor(pca1_scores[common_score_cols],
                                    med_scores[common_score_cols], method = 'spearman')
      } else {
        rho_mean_pca1 <- 0
        rho_pca1_med <- 0
      }

      radar_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_med'] <- rho_mean_med
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['rho_pca1_med'] <- rho_pca1_med
      radar_values[[names_sigs[k]]][[names_datasets[i]]]['rho_mean_pca1'] <- rho_mean_pca1

      if (length(pca1_scores) > 1) {
        bars_plot <- props_of_variances[1:min(10, length(props_of_variances))]
        radar_values[[names_sigs[k]]][[names_datasets[i]]]['prop_pca1_var'] <- bars_plot[1]
      } else {
        radar_values[[names_sigs[k]]][[names_datasets[i]]]['prop_pca1_var'] <- 0
      }

      # --- Build combined score correlation matrix ---
      # Expand common_score_cols to include ES scores
      if (!is.null(es.gsva)) {
        es_common <- intersect(names(es.gsva), intersect(names(es.ssGSEA), names(es.plage)))
        if (!is.null(pca1_scores)) {
          all_common <- intersect(es_common,
                                   intersect(common_score_cols, names(pca1_scores)))
          output_mat <- cbind(mean_scores[all_common], med_scores[all_common],
                               pca1_scores[all_common], es.ssGSEA[all_common],
                               es.gsva[all_common], es.plage[all_common])
          colnames(output_mat) <- c('Mean_Scores', 'Median_Scores', 'PCA1_Scores',
                                    'ssGSEA', 'GSVA', 'PLAGE')
        } else {
          all_common <- intersect(es_common, common_score_cols)
          output_mat <- cbind(mean_scores[all_common], med_scores[all_common],
                               es.ssGSEA[all_common], es.gsva[all_common], es.plage[all_common])
          colnames(output_mat) <- c('Mean_Scores', 'Median_Scores', 'ssGSEA', 'GSVA', 'PLAGE')
        }
        if (length(all_common) > 1) {
          score_cor_mats[[paste0(names_datasets[i], '_', names_sigs[k])]] <-
            stats::cor(output_mat, method = 'spearman')
        }
      }

      # --- Mixture models ---
      mixture_models[[names_sigs[k]]][[names_datasets[i]]] <- list(
        median = NULL, mean = NULL, pca1 = NULL
      )

      if (!radar_only) {
        if (length(med_scores) > 1) {
          max_clusters <- min(ceiling(sum(!is.na(med_scores)) / 2), 10)
          mixture_models[[names_sigs[k]]][[names_datasets[i]]][['median']] <-
            mclust::Mclust(med_scores, G = 1:max_clusters)
        }
        if (length(mean_scores) > 1) {
          max_clusters <- min(ceiling(sum(!is.na(mean_scores)) / 2), 10)
          mixture_models[[names_sigs[k]]][[names_datasets[i]]][['mean']] <-
            mclust::Mclust(mean_scores, G = 1:max_clusters)
        }
        if (length(pca1_scores) > 1) {
          max_clusters <- min(ceiling(sum(!is.na(pca1_scores)) / 2), 10)
          mixture_models[[names_sigs[k]]][[names_datasets[i]]][['pca1']] <-
            mclust::Mclust(pca1_scores, G = 1:max_clusters)
        }
      }

      # Store all scores for plotting
      scores_all[[names_sigs[k]]][[names_datasets[i]]] <- list(
        med_scores = med_scores,
        mean_scores = mean_scores,
        pca1_scores = pca1_scores,
        es.gsva = es.gsva,
        es.ssGSEA = es.ssGSEA,
        es.plage = es.plage,
        common_score_cols = common_score_cols,
        props_of_variances = props_of_variances
      )
    }
  }

  list(
    radar_values = radar_values,
    scores = scores_all,
    pca_results = pca_results,
    score_cor_mats = score_cor_mats,
    mixture_models = mixture_models
  )
}

# ============================================================================
# Helper: format mixture model output string
# ============================================================================
.format_mixture_string <- function(model, sig_name, dataset_name, metric_name) {
  if (is.null(model)) return(NULL)
  output_string <- paste0(sig_name, ' ', dataset_name, ', ', metric_name, ' score: There ')
  if (model$G == 1) {
    output_string <- paste0(output_string, ' is one component in the Gaussian mixture model. ')
  } else {
    output_string <- paste0(output_string, ' are ', model$G,
                            ' components in the Gaussian mixture model. ')
  }
  if (model$modelName == 'E') {
    output_string <- paste0(output_string,
                            'Best model is Gaussian distributions with equal variance (E). ')
  } else if (model$modelName == 'V') {
    output_string <- paste0(output_string,
                            'Best model is Gaussian distributions with variable variances (V). ')
  } else {
    output_string <- paste0(output_string, 'Best model is univariate Gaussian distribution. ')
  }
  paste0(output_string, '\n')
}

# ============================================================================
# plot_metrics: Side effects — PDFs + CSV files
# ============================================================================
plot_metrics <- function(compute_result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
                         names_datasets, out_dir, file = NULL, showResults = FALSE) {

  if (!dir.exists(file.path(out_dir, 'metrics_tables'))) {
    dir.create(file.path(out_dir, 'metrics_tables'))
  }

  # ---- Section 1: Mean/Median/PCA1 scatter plots ----
  for (k in seq_along(names_sigs)) {
    max_title_length <- max(nchar(paste0(names_datasets, ' ', names_sigs[k])))

    if (showResults) grDevices::dev.new()
    else grDevices::pdf(file.path(out_dir, paste0('sig_compare_metrics_', names_sigs[k], '.pdf')),
                        width = 3 * length(names_datasets), height = 10)
    graphics::par(mfcol = c(4, length(names_datasets)), mar = c(4, 4, 4, 4))

    for (i in seq_along(names_datasets)) {
      sc <- compute_result$scores[[names_sigs[k]]][[names_datasets[i]]]
      med_scores <- sc$med_scores
      mean_scores <- sc$mean_scores
      pca1_scores <- sc$pca1_scores
      common_score_cols <- sc$common_score_cols
      props_of_variances <- sc$props_of_variances

      jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

      # Plot 1: Median vs Mean
      if (length(common_score_cols) > 1) {
        graphics::smoothScatter(med_scores[common_score_cols], mean_scores[common_score_cols],
                                colramp = jet.colors, xlab = NA, ylab = NA, main = 'Median vs Mean')
        graphics::points(med_scores[common_score_cols], mean_scores[common_score_cols], pch = '.')
        graphics::par(new = TRUE)
        graphics::plot(stats::density(med_scores[common_score_cols]), axes = FALSE,
                       xlab = NA, ylab = NA, col = 'red', main = NA, lwd = 2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
        graphics::mtext(side = 2, line = 2, 'Mean', cex = 0.8)
        graphics::mtext(side = 1, line = 2, 'Median', cex = 0.8)
        graphics::mtext(side = 3, line = 2.5, paste0(names_datasets[i], ' ', names_sigs[k]),
                        cex = min(1, 3 * 10 / max_title_length))
        rho <- stats::cor(med_scores[common_score_cols], mean_scores[common_score_cols], method = 'spearman')
        graphics::mtext(paste0('rho = ', format(rho, digits = 2)), side = 3, line = 0,
                        cex = 0.6, at = max(med_scores[common_score_cols]))
      } else {
        graphics::plot.new()
        graphics::mtext(side = 3, line = 2.5, paste0(names_datasets[i], ' ', names_sigs[k]),
                        cex = min(1, 3 * 10 / max_title_length))
        graphics::title(paste0('\n\nToo many NA values for Mean/Median in \n',
                               names_datasets[i], ' ', names_sigs[k]))
      }

      # Plot 2: Mean vs PCA1
      if (length(pca1_scores) > 1) {
        graphics::smoothScatter(mean_scores[common_score_cols], pca1_scores[common_score_cols],
                                colramp = jet.colors, xlab = NA, ylab = NA, main = 'Mean vs PCA1')
        graphics::points(mean_scores[common_score_cols], pca1_scores[common_score_cols], pch = '.')
        graphics::par(new = TRUE)
        graphics::plot(stats::density(mean_scores[common_score_cols]), axes = FALSE,
                       xlab = NA, ylab = NA, col = 'red', main = NA, lwd = 2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
        rho <- stats::cor(mean_scores[common_score_cols], pca1_scores[common_score_cols], method = 'spearman')
        graphics::mtext(paste0('rho = ', format(rho, digits = 2)), side = 3, line = 0,
                        cex = 0.6, at = max(mean_scores[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'PCA1', cex = 0.8)
        graphics::mtext(side = 1, line = 2, 'Mean', cex = 0.8)

        # Plot 3: PCA1 vs Median
        graphics::smoothScatter(pca1_scores[common_score_cols], med_scores[common_score_cols],
                                colramp = jet.colors, xlab = NA, ylab = NA, main = 'PCA1 vs Median')
        graphics::points(pca1_scores[common_score_cols], med_scores[common_score_cols], pch = '.')
        graphics::par(new = TRUE)
        graphics::plot(stats::density(pca1_scores[common_score_cols]), axes = FALSE,
                       xlab = NA, ylab = NA, col = 'red', main = NA, lwd = 2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
        rho <- stats::cor(pca1_scores[common_score_cols], med_scores[common_score_cols], method = 'spearman')
        graphics::mtext(paste0('rho = ', format(rho, digits = 2)), side = 3, line = 0,
                        cex = 0.6, at = max(pca1_scores[common_score_cols]))
        graphics::mtext(side = 2, line = 2, 'Median', cex = 0.8)
        graphics::mtext(side = 1, line = 2, 'PCA1', cex = 0.8)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1/Mean in \n',
                               names_datasets[i], ' ', names_sigs[k]))
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Median/PCA1 in \n',
                               names_datasets[i], ' ', names_sigs[k]))
      }

      # Plot 4: Scree plot
      if (length(pca1_scores) > 1) {
        bars_plot <- props_of_variances[1:min(10, length(props_of_variances))]
        output_mat <- as.matrix(props_of_variances)
        colnames(output_mat) <- c('Proportion of variance')
        row.names(output_mat) <- paste0('PCA ', seq_len(nrow(output_mat)))
        utils::write.table(output_mat,
                           file = file.path(out_dir, 'metrics_tables',
                                            paste0('pca_vs_var_', names_sigs[k], '_', names_datasets[i], '.txt')),
                           quote = FALSE, sep = '\t')
        graphics::barplot(bars_plot, main = "PCA vs proportion\n of variance")
        graphics::mtext(side = 1, line = 2, 'PCA', cex = 0.8)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1 in \n',
                               names_datasets[i], ' ', names_sigs[k]))
      }

      # Write PCA loadings
      pca_obj <- compute_result$pca_results[[names_sigs[k]]][[names_datasets[i]]]$pca_obj
      if (!is.null(pca_obj)) {
        output_mat <- pca_obj$x
        data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
        if (!is.null(colnames(data.matrix))) row.names(output_mat) <- colnames(data.matrix)
        utils::write.table(output_mat,
                           file = file.path(out_dir, 'metrics_tables',
                                            paste0('pca_loadings_', names_sigs[k], '_', names_datasets[i], '.txt')),
                           quote = FALSE, sep = '\t')
      }
    }

    if (showResults) {
      grDevices::dev.copy(grDevices::pdf, file.path(out_dir, paste0('sig_compare_metrics_', names_sigs[k], '.pdf')),
                          width = 3 * length(names_datasets), height = 10)
    }
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
  }
  cat('Metrics compared successfully.\n', file = file)

  # ---- Section 2: ES scatter plots ----
  for (k in seq_along(names_sigs)) {
    max_title_length <- max(nchar(paste0(names_datasets, ' ', names_sigs[k])))

    if (showResults) grDevices::dev.new()
    else grDevices::pdf(file.path(out_dir, paste0('sig_compare_ES_metrics_', names_sigs[k], '.pdf')),
                        width = 3 * length(names_datasets), height = 7.5)
    graphics::par(mfcol = c(3, length(names_datasets)), mar = c(4, 4, 4, 4))

    for (i in seq_along(names_datasets)) {
      sc <- compute_result$scores[[names_sigs[k]]][[names_datasets[i]]]
      es.gsva <- sc$es.gsva
      es.ssGSEA <- sc$es.ssGSEA
      es.plage <- sc$es.plage

      if (!is.null(es.gsva)) {
        common_es <- intersect(names(es.gsva), intersect(names(es.ssGSEA), names(es.plage)))
      } else {
        common_es <- character(0)
      }

      jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

      # GSVA vs ssGSEA
      if (length(common_es) > 1) {
        graphics::smoothScatter(es.gsva[common_es], es.ssGSEA[common_es],
                                colramp = jet.colors, xlab = NA, ylab = NA, main = 'GSVA vs ssGSEA')
        graphics::points(es.gsva[common_es], es.ssGSEA[common_es], pch = '.')
        graphics::par(new = TRUE)
        graphics::plot(stats::density(es.gsva[common_es]), axes = FALSE,
                       xlab = NA, ylab = NA, col = 'red', main = NA, lwd = 2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
        graphics::mtext(side = 2, line = 2, 'ssGSEA', cex = 0.8)
        graphics::mtext(side = 1, line = 2, 'GSVA', cex = 0.8)
        graphics::mtext(side = 3, line = 2.5, paste0(names_datasets[i], ' ', names_sigs[k]),
                        cex = min(1, 3 * 10 / max_title_length))
        rho <- stats::cor(es.gsva[common_es], es.ssGSEA[common_es], method = 'spearman')
        graphics::mtext(paste0('rho = ', format(rho, digits = 2)), side = 3, line = 0,
                        cex = 0.6, at = max(es.gsva[common_es]))
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for ssGSEA/GSVA'))
      }

      # ssGSEA vs PLAGE
      if (length(common_es) > 1) {
        graphics::smoothScatter(es.ssGSEA[common_es], es.plage[common_es],
                                colramp = jet.colors, xlab = NA, ylab = NA, main = 'ssGSEA vs PLAGE')
        graphics::points(es.ssGSEA[common_es], es.plage[common_es], pch = '.')
        graphics::par(new = TRUE)
        graphics::plot(stats::density(es.ssGSEA[common_es]), axes = FALSE,
                       xlab = NA, ylab = NA, col = 'red', main = NA, lwd = 2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
        rho <- stats::cor(es.ssGSEA[common_es], es.plage[common_es], method = 'spearman')
        graphics::mtext(paste0('rho = ', format(rho, digits = 2)), side = 3, line = 0,
                        cex = 0.6, at = max(es.ssGSEA[common_es]))
        graphics::mtext(side = 2, line = 2, 'PLAGE', cex = 0.8)
        graphics::mtext(side = 1, line = 2, 'ssGSEA', cex = 0.8)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PLAGE/ssGSEA'))
      }

      # PLAGE vs GSVA
      if (length(common_es) > 1) {
        graphics::smoothScatter(es.plage[common_es], es.gsva[common_es],
                                colramp = jet.colors, xlab = NA, ylab = NA, main = 'PLAGE vs GSVA')
        graphics::points(es.plage[common_es], es.gsva[common_es], pch = '.')
        graphics::par(new = TRUE)
        graphics::plot(stats::density(es.plage[common_es]), axes = FALSE,
                       xlab = NA, ylab = NA, col = 'red', main = NA, lwd = 2)
        graphics::axis(side = 4)
        graphics::mtext(side = 4, line = 2, 'Density', cex = 0.8)
        rho <- stats::cor(es.plage[common_es], es.gsva[common_es], method = 'spearman')
        graphics::mtext(paste0('rho = ', format(rho, digits = 2)), side = 3, line = 0,
                        cex = 0.6, at = max(es.plage[common_es]))
        graphics::mtext(side = 2, line = 2, 'GSVA', cex = 0.8)
        graphics::mtext(side = 1, line = 2, 'PLAGE', cex = 0.8)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for GSVA/PLAGE'))
      }

      # Write combined metrics table
      sc2 <- compute_result$scores[[names_sigs[k]]][[names_datasets[i]]]
      pca1 <- sc2$pca1_scores
      csc <- sc2$common_score_cols
      if (!is.null(es.gsva) && length(common_es) > 1) {
        if (!is.null(pca1)) {
          all_c <- intersect(common_es, intersect(csc, names(pca1)))
          if (length(all_c) > 1) {
            out_mat <- cbind(sc2$mean_scores[all_c], sc2$med_scores[all_c], pca1[all_c],
                             es.ssGSEA[all_c], es.gsva[all_c], es.plage[all_c])
            colnames(out_mat) <- c('Mean_Scores', 'Median_Scores', 'PCA1_Scores',
                                    'ssGSEA', 'GSVA', 'PLAGE')
            utils::write.table(out_mat,
                               file = file.path(out_dir, 'metrics_tables',
                                                paste0('metrics_table_', names_sigs[k], '_', names_datasets[i], '.txt')),
                               quote = FALSE, sep = '\t')
          }
        } else {
          all_c <- intersect(common_es, csc)
          if (length(all_c) > 1) {
            out_mat <- cbind(sc2$mean_scores[all_c], sc2$med_scores[all_c],
                             es.ssGSEA[all_c], es.gsva[all_c], es.plage[all_c])
            colnames(out_mat) <- c('Mean_Scores', 'Median_Scores', 'ssGSEA', 'GSVA', 'PLAGE')
            utils::write.table(out_mat,
                               file = file.path(out_dir, 'metrics_tables',
                                                paste0('metrics_table_', names_sigs[k], '_', names_datasets[i], '.txt')),
                               quote = FALSE, sep = '\t')
          }
        }
      }
    }

    if (showResults) {
      grDevices::dev.copy(grDevices::pdf, file.path(out_dir, paste0('sig_compare_ES_metrics_', names_sigs[k], '.pdf')),
                          width = 3 * length(names_datasets), height = 7.5)
    }
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
  }
  cat('ES scores computed successfully.\n', file = file)

  # ---- Section 3: Scoring metric correlation heatmaps ----
  nice_row_names <- c('Mean', 'Median', 'PCA1', 'GSVA', 'ssGSEA', 'PLAGE')
  names(nice_row_names) <- c('Mean_Scores', 'Median_Scores', 'PCA1_Scores', 'GSVA', 'ssGSEA', 'PLAGE')

  for (k in seq_along(names_sigs)) {
    for (i in seq_along(names_datasets)) {
      key <- paste0(names_datasets[i], '_', names_sigs[k])
      plot_mat <- compute_result$score_cor_mats[[key]]
      if (is.null(plot_mat)) next

      if (showResults) grDevices::dev.new()
      else grDevices::pdf(file.path(out_dir, paste0('scoring_metrics_corr_', names_datasets[i], '_', names_sigs[k], '.pdf')),
                          width = 4, height = 4)
      graphics::par(mfcol = c(4, length(names_datasets)), mar = c(4, 4, 4, 4))

      row.names(plot_mat) <- nice_row_names[rownames(plot_mat)]
      ans_hmap <- ComplexHeatmap::Heatmap(plot_mat, show_column_dend = FALSE,
                                           show_column_names = FALSE,
                                           name = names_datasets[i],
                                           col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                                           heatmap_legend_param = list(title = 'Correlation', color_bar = "continuous",
                                                                       legend_direction = 'vertical'),
                                           column_title = paste0(names_datasets[i], ' ', names_sigs[k],
                                                                 '\nScoring Metric Correlation'),
                                           row_names_gp = grid::gpar(fontsize = 10),
                                           row_title = 'Scoring metrics')
      ComplexHeatmap::draw(ans_hmap, heatmap_legend_side = "left")

      if (showResults) {
        grDevices::dev.copy(grDevices::pdf, file.path(out_dir, paste0('scoring_metrics_corr_', names_datasets[i], '_', names_sigs[k], '.pdf')),
                            width = 7, height = 3.5)
      }
      if (grDevices::dev.cur() != 1) grDevices::dev.off()
    }
  }

  # ---- Section 4: Gaussian mixture model BIC plots ----
  mixture_model.out <- file.path(out_dir, "mixture_model_out.txt")
  mixture_model.con <- file(mixture_model.out, open = "w")

  for (k in seq_along(names_sigs)) {
    max_title_length <- max(nchar(paste0(names_datasets, ' ', names_sigs[k])))

    if (showResults) grDevices::dev.new()
    else grDevices::pdf(file.path(out_dir, paste0('sig_gaussian_mixture_model_', names_sigs[k], '.pdf')),
                        width = 3 * length(names_datasets), height = 10)
    graphics::par(mfcol = c(3, length(names_datasets)), mar = c(4, 4, 4, 4))

    for (i in seq_along(names_datasets)) {
      mm <- compute_result$mixture_models[[names_sigs[k]]][[names_datasets[i]]]

      # Median BIC
      if (!is.null(mm[['median']])) {
        mclust::plot.Mclust(x = mm[['median']], what = 'BIC', main = 'Median score')
        graphics::mtext(side = 3, line = 2.5, paste0(names_datasets[i], ' ', names_sigs[k]),
                        cex = min(1, 3 * 10 / max_title_length))
        cat(.format_mixture_string(mm[['median']], names_sigs[k], names_datasets[i], 'Median'),
            file = mixture_model.con)
      } else {
        graphics::plot.new()
        graphics::mtext(side = 3, line = 2.5, paste0(names_datasets[i], ' ', names_sigs[k]),
                        cex = min(1, 3 * 10 / max_title_length))
        graphics::title(paste0('\n\nToo many NA values for Median in \n',
                               names_datasets[i], ' ', names_sigs[k]))
      }

      # Mean BIC
      if (!is.null(mm[['mean']])) {
        mclust::plot.Mclust(x = mm[['mean']], what = 'BIC', main = 'Mean score')
        cat(.format_mixture_string(mm[['mean']], names_sigs[k], names_datasets[i], 'Mean'),
            file = mixture_model.con)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Mean in \n',
                               names_datasets[i], ' ', names_sigs[k]))
      }

      # PCA1 BIC
      if (!is.null(mm[['pca1']])) {
        mclust::plot.Mclust(x = mm[['pca1']], what = 'BIC', main = 'PCA1 score')
        cat(.format_mixture_string(mm[['pca1']], names_sigs[k], names_datasets[i], 'PCA1'),
            file = mixture_model.con)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1 in \n',
                               names_datasets[i], ' ', names_sigs[k]))
      }
    }

    if (showResults) {
      grDevices::dev.copy(grDevices::pdf, file.path(out_dir, paste0('sig_gaussian_mixture_model_', names_sigs[k], '.pdf')),
                          width = 3 * length(names_datasets), height = 10)
    }
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
  }
  close(mixture_model.con)
  save(list = "mixture_models", envir = list2env(list(mixture_models = compute_result$mixture_models)),
       file = file.path(out_dir, "mixture_models_raw_out.rda"))
  cat('Gaussian mixture models computed successfully.\n', file = file)

  # ---- Section 5: QQ plots ----
  for (k in seq_along(names_sigs)) {
    max_title_length <- max(nchar(paste0(names_datasets, ' ', names_sigs[k])))

    if (showResults) grDevices::dev.new()
    else grDevices::pdf(file.path(out_dir, paste0('sig_qq_plots_', names_sigs[k], '.pdf')),
                        width = 3 * length(names_datasets), height = 10)
    graphics::par(mfcol = c(3, length(names_datasets)), mar = c(4, 4, 4, 4))

    for (i in seq_along(names_datasets)) {
      sc <- compute_result$scores[[names_sigs[k]]][[names_datasets[i]]]

      if (length(sc$med_scores) > 1) {
        stats::qqnorm(sc$med_scores, plot.it = TRUE, main = 'Median score')
        graphics::mtext(side = 3, line = 2.5, paste0(names_datasets[i], ' ', names_sigs[k]),
                        cex = min(1, 3 * 10 / max_title_length))
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Median'))
      }

      if (length(sc$mean_scores) > 1) {
        stats::qqnorm(sc$mean_scores, plot.it = TRUE, main = 'Mean score')
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for Mean'))
      }

      if (length(sc$pca1_scores) > 1) {
        stats::qqnorm(sc$pca1_scores, plot.it = TRUE, main = 'PCA1 score')
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\nToo many NA values for PCA1'))
      }
    }

    if (showResults) {
      grDevices::dev.copy(grDevices::pdf, file.path(out_dir, paste0('sig_qq_plots_', names_sigs[k], '.pdf')),
                          width = 3 * length(names_datasets), height = 10)
    }
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
  }
  cat('QQ plots computed successfully.\n', file = file)
}

# ============================================================================
# compare_metrics_loc: Original API wrapper — backward-compatible
# ============================================================================
compare_metrics_loc <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                                 out_dir = '~', file = NULL, showResults = FALSE, radar_plot_values) {

  result <- compute_metrics(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
  plot_metrics(result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
               names_datasets, out_dir, file, showResults)

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
