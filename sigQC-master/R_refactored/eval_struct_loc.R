# eval_struct_loc.R — REFACTORED
#
# Evaluates gene signature structure via hierarchical clustering and biclustering.
# Split into compute_struct() + plot_struct(). The duplicated z-transform/binarize/biclust
# pipeline (repeated 4x in the original) is now computed once in compute_struct().

# ============================================================================
# Helper: z-transform + binarize + bicluster for one sig-dataset pair
# ============================================================================
.compute_biclust <- function(data.matrix, inter) {
  sig_scores <- as.matrix(data.matrix[inter, ])
  sig_scores[!is.finite(sig_scores)] <- NA

  # Z-transform each gene (with zero-variance guard)
  for (gene in inter) {
    gene_values <- as.numeric(sig_scores[gene, ])
    gene_sd <- stats::sd(gene_values, na.rm = TRUE)
    if (is.na(gene_sd) || gene_sd == 0) {
      sig_scores[gene, ] <- 0
    } else {
      sig_scores[gene, ] <- (gene_values - mean(gene_values, na.rm = TRUE)) / gene_sd
    }
  }

  # Binarize at midpoint of z-score range
  z_vals <- stats::na.omit(t(sig_scores))
  threshold <- min(z_vals) + (max(z_vals) - min(z_vals)) / 2
  x <- z_vals
  x[x <= threshold] <- 0
  x[x > threshold] <- 1

  # Run BCCC biclustering
  Xmotif <- biclust::biclust(x, method = biclust::BCCC(), delta = 1, alpha = 1.5, number = 50)

  list(
    z_scores = sig_scores,
    binarized = x,
    biclust_result = Xmotif,
    threshold = threshold
  )
}

# ============================================================================
# compute_struct: Pure computation — no side effects
# ============================================================================
# Returns a list with:
#   $sig_scores_all_mats — padded expression matrices per sig-dataset
#   $all_row_names       — union of gene names per signature
#   $biclust_results     — nested [[sig]][[dataset]] with z_scores, binarized, biclust_result
#   $any_biclusters      — logical, whether any sig-dataset pair has >1 bicluster
compute_struct <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                           covariates = NULL) {
  # --- Step 1: Build union of gene names per signature ---
  all_row_names <- list()
  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    all_row_names[[names_sigs[k]]] <- c()
    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      all_row_names[[names_sigs[k]]] <- c(all_row_names[[names_sigs[k]]], inter)
    }
    all_row_names[[names_sigs[k]]] <- unique(all_row_names[[names_sigs[k]]])
  }

  # --- Step 2: Build padded expression matrices ---
  sig_scores_all_mats <- list()
  for (k in seq_along(names_sigs)) {
    sig_scores_all_mats[[names_sigs[k]]] <- list()
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      sig_scores <- as.matrix(data.matrix[inter, ])
      rows_needed <- setdiff(all_row_names[[names_sigs[k]]], inter)
      if (length(rows_needed) > 0) {
        sig_scores <- rbind(sig_scores, matrix(NA, nrow = length(rows_needed), ncol = ncol(sig_scores)))
        row.names(sig_scores) <- c(inter, rows_needed)
      }
      sig_scores_all_mats[[names_sigs[k]]][[names_datasets[i]]] <- sig_scores
    }
  }

  # --- Step 3: Run biclustering once per sig-dataset pair ---
  biclust_results <- list()
  any_biclusters <- FALSE
  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    biclust_results[[names_sigs[k]]] <- list()
    for (i in seq_along(names_datasets)) {
      data.matrix <- mRNA_expr_matrix[[names_datasets[i]]]
      inter <- intersect(gene_sig, row.names(data.matrix))
      bc <- .compute_biclust(data.matrix, inter)
      biclust_results[[names_sigs[k]]][[names_datasets[i]]] <- bc
      if (bc$biclust_result@Number > 1) any_biclusters <- TRUE
    }
  }

  list(
    sig_scores_all_mats = sig_scores_all_mats,
    all_row_names = all_row_names,
    biclust_results = biclust_results,
    any_biclusters = any_biclusters
  )
}

# ============================================================================
# plot_struct: Side effects — PDFs
# ============================================================================
plot_struct <- function(compute_result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
                        names_datasets, covariates, out_dir, file = NULL, showResults = FALSE) {
  max_title_length <- max(nchar(outer(names_datasets, names_sigs, paste)))

  # ---- Section 1: Hierarchical clustering heatmaps ----
  for (k in seq_along(names_sigs)) {
    gene_sig <- gene_sigs_list[[names_sigs[k]]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)

    hmaps <- sapply(seq_along(names_datasets), function(i) {
      sig_scores <- compute_result$sig_scores_all_mats[[names_sigs[k]]][[names_datasets[i]]]
      tryCatch({
        dim.pdf <- dim(sig_scores)
        h <- dim.pdf[1]
        row_names.fontsize <- if (h < 20) 12 else 5 / log10(h)

        if (length(covariates[[names_datasets[i]]]) == 0) {
          ans_hmap <- ComplexHeatmap::Heatmap(sig_scores, show_column_dend = FALSE,
                                               show_column_names = FALSE, name = names_datasets[i],
                                               heatmap_legend_param = list(title = names_datasets[i],
                                                                           color_bar = "continuous",
                                                                           legend_direction = 'vertical'),
                                               column_title = names_datasets[i],
                                               row_names_gp = grid::gpar(fontsize = row_names.fontsize),
                                               row_title = 'Genes')
        } else {
          if (is.vector(covariates[[names_datasets[i]]][['annotations']])) {
            ha1 <- ComplexHeatmap::HeatmapAnnotation(
              df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][
                intersect(names(covariates[[names_datasets[i]]][['annotations']]), colnames(sig_scores))]),
              col = covariates[[names_datasets[i]]][['colors']], na_col = "grey")
          } else {
            ha1 <- ComplexHeatmap::HeatmapAnnotation(
              df = as.data.frame(covariates[[names_datasets[i]]][['annotations']][
                intersect(rownames(covariates[[names_datasets[i]]][['annotations']]), colnames(sig_scores)), ]),
              col = covariates[[names_datasets[i]]][['colors']], na_col = "grey", which = "column")
          }
          ans_hmap <- ComplexHeatmap::Heatmap(sig_scores, show_column_dend = FALSE,
                                               show_column_names = FALSE, name = names_datasets[i],
                                               heatmap_legend_param = list(title = names_datasets[i],
                                                                           color_bar = "continuous",
                                                                           legend_direction = 'vertical'),
                                               column_title = names_datasets[i],
                                               row_names_gp = grid::gpar(fontsize = row_names.fontsize),
                                               row_title = 'Genes', top_annotation = ha1)
        }
        ans_hmap
      }, error = function(err) {
        cat(paste0('Error when creating expression heatmaps for ', names_datasets[i], ' ',
                   names_sigs[k], ': ', err, '\n'), file = file)
        graphics::plot.new()
        graphics::title(paste0('\n \n \n', names_datasets[i], ' ', names_sigs[k]))
        NULL
      })
    })

    for (i in seq_along(names_datasets)) {
      all_hmaps <- hmaps[[1]]
      if (length(hmaps) > 1) {
        lapply(hmaps[2:length(hmaps)], function(x) all_hmaps <<- ComplexHeatmap::add_heatmap(all_hmaps, x))
      }
      if (showResults) grDevices::dev.new()
      else grDevices::pdf(file.path(out_dir, paste0('sig_eval_struct_clustering_', names_datasets[i], '_', names_sigs[k], '.pdf')),
                          width = 5 * length(names_datasets), height = 10)
      ComplexHeatmap::draw(all_hmaps, heatmap_legend_side = "left",
                           annotation_legend_side = "left", main_heatmap = names_datasets[i])
      if (showResults) {
        grDevices::dev.copy(grDevices::pdf, file.path(out_dir, paste0('sig_eval_struct_clustering_', names_datasets[i], '_', names_sigs[k], '.pdf')),
                            width = 5 * length(names_datasets), height = 10)
      }
      cat('Expression heatmaps saved successfully.\n', file = file)
      if (grDevices::dev.cur() != 1) grDevices::dev.off()
    }
  }

  # ---- Section 2: Biclustering heatmaps (continuous z-scores) ----
  if (compute_result$any_biclusters) {
    grDevices::dev.new()
    graphics::par(cex.main = 0.8, cex.lab = 0.8, oma = c(4, 2, 2, 2), mar = c(4, 4, 4, 4))

    hmaps <- lapply(seq_len(length(names_datasets) * length(names_sigs)), function(idx) {
      dataset_ind <- ((idx - 1) %% length(names_datasets)) + 1
      sig_ind <- ceiling(idx / length(names_datasets))
      bc <- compute_result$biclust_results[[names_sigs[sig_ind]]][[names_datasets[dataset_ind]]]

      if (bc$biclust_result@Number > 1) {
        biclust::heatmapBC(stats::na.omit(t(bc$z_scores)), bicResult = bc$biclust_result,
                           col = gplots::colorpanel(100, "blue", "white", "red"),
                           xlab = 'Gene ID', ylab = 'Sample')
        graphics::title(paste0('\n \n \nBivariate clustering\n',
                               names_datasets[dataset_ind], ' ', names_sigs[sig_ind]),
                        cex = min(1, 4 * 10 / max_title_length))
        graphics::axis(1, at = seq_along(rownames(bc$z_scores)),
                       labels = rownames(bc$z_scores), las = 2, tck = 0, cex.axis = 0.6)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\n\n <=1 bivariate clusters for\n',
                               names_datasets[dataset_ind], ' ', names_sigs[sig_ind]),
                        cex = min(1, 4 * 10 / max_title_length))
        cat(paste0('<= 1 bi-clusters found for: ', names_datasets[dataset_ind], ' ',
                   names_sigs[sig_ind], '\n'), file = file)
      }
      grab_grob()
    })

    draw.heatmaps(hmaps, names_datasets, names_sigs)
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_eval_bivariate_clustering.pdf'),
                        width = 4 * length(names_datasets), height = 4 * length(names_sigs))
    cat('Bi-clustering completed successfully\n', file = file)
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
  } else {
    cat('Bi-clustering completed successfully. No bi-clusters found.\n', file = file)
  }

  # ---- Section 3: Biclustering heatmaps (binarized 0/1) ----
  if (compute_result$any_biclusters) {
    grDevices::dev.new()
    graphics::par(cex.main = 0.8, cex.lab = 0.8, oma = c(4, 2, 2, 2), mar = c(4, 4, 4, 4))

    hmaps <- lapply(seq_len(length(names_datasets) * length(names_sigs)), function(idx) {
      dataset_ind <- ((idx - 1) %% length(names_datasets)) + 1
      sig_ind <- ceiling(idx / length(names_datasets))
      bc <- compute_result$biclust_results[[names_sigs[sig_ind]]][[names_datasets[dataset_ind]]]

      if (bc$biclust_result@Number > 1) {
        biclust::heatmapBC(bc$binarized, bicResult = bc$biclust_result,
                           col = gplots::colorpanel(100, "blue", "white", "red"),
                           xlab = 'Gene ID', ylab = 'Sample')
        graphics::title(paste0('\n \n \nBivariate clustering\n',
                               names_datasets[dataset_ind], ' ', names_sigs[sig_ind]),
                        cex = min(1, 4 * 10 / max_title_length))
        graphics::axis(1, at = seq_along(rownames(bc$z_scores)),
                       labels = rownames(bc$z_scores), las = 2, tck = 0, cex.axis = 0.6)
      } else {
        graphics::plot.new()
        graphics::title(paste0('\n\n\n <=1 bivariate clusters for\n',
                               names_datasets[dataset_ind], ' ', names_sigs[sig_ind]),
                        cex = min(1, 4 * 10 / max_title_length))
      }
      grab_grob()
    })

    draw.heatmaps(hmaps, names_datasets, names_sigs)
    grDevices::dev.copy(grDevices::pdf, file.path(out_dir, 'sig_eval_bivariate_clustering_binarized_maps.pdf'),
                        width = 4 * length(names_datasets), height = 4 * length(names_sigs))
    if (grDevices::dev.cur() != 1) grDevices::dev.off()
  } else {
    cat('Binarized bi-clustering completed. No binarized bi-clusters found.\n', file = file)
  }
  if (grDevices::dev.cur() != 1) grDevices::dev.off()
}

# ============================================================================
# eval_struct_loc: Original API wrapper — backward-compatible
# ============================================================================
eval_struct_loc <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets,
                            covariates, out_dir = '~', file = NULL, showResults = FALSE,
                            radar_plot_values) {

  result <- compute_struct(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets, covariates)
  plot_struct(result, gene_sigs_list, names_sigs, mRNA_expr_matrix,
              names_datasets, covariates, out_dir, file, showResults)

  # eval_struct_loc doesn't add any radar chart metrics
  radar_plot_values
}
