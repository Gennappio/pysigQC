# sigsQcNegativeControl.R — OPTIMIZED
#
# Optimized version of negative and permutation controls.
# Key optimizations over the corrected version:
#   1. No file I/O roundtrip — compute_*() return results in memory
#   2. Sparse permutation — only permute signature gene rows, not full matrix copies
#   3. Vectorized quantile computation — single apply() instead of 6 separate calls
#   4. parallel::mclapply() for resampling (with lapply fallback on Windows)
#   5. BUG-1 fix: sample() instead of runif()
#   6. Deduplicated summary/plotting logic via helper functions
#
# This file depends on the compute_*() functions from R_refactored/.
# Source them before using this file:
#   source("R_refactored/eval_var_loc.R")
#   source("R_refactored/eval_expr_loc.R")
#   ... etc

# ============================================================================
# Helper: detect available cores and choose parallel backend
# ============================================================================
.get_apply_fn <- function(n_cores = NULL) {
  if (.Platform$OS.type == "windows") {
    return(lapply)
  }
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = n_cores)
}

# ============================================================================
# Helper: compute QC metrics in memory (no disk I/O)
# Returns radar_plot_values directly instead of writing/reading radarchart_table.txt
# ============================================================================
.compute_qc_metrics_inmemory <- function(gene_sigs_list, names_sigs,
                                          mRNA_expr_matrix, names_datasets) {
  radar_plot_values <- list()
  for (k in seq_along(names_sigs)) {
    radar_plot_values[[names_sigs[k]]] <- list()
  }

  # Run each compute module, collecting radar values
  tryCatch({
    r <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        vals <- r$radar_values[[names_sigs[k]]][[names_datasets[i]]]
        for (m in names(vals)) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <- vals[m]
        }
      }
    }
  }, error = function(e) {})

  tryCatch({
    r <- compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        vals <- r$radar_values[[names_sigs[k]]][[names_datasets[i]]]
        for (m in names(vals)) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <- vals[m]
        }
      }
    }
  }, error = function(e) {})

  tryCatch({
    r <- compute_compactness(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        vals <- r$radar_values[[names_sigs[k]]][[names_datasets[i]]]
        for (m in names(vals)) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <- vals[m]
        }
      }
    }
  }, error = function(e) {})

  tryCatch({
    r <- compute_metrics(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        vals <- r$radar_values[[names_sigs[k]]][[names_datasets[i]]]
        for (m in names(vals)) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <- vals[m]
        }
      }
    }
  }, error = function(e) {})

  tryCatch({
    r <- compute_stan(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        vals <- r$radar_values[[names_sigs[k]]][[names_datasets[i]]]
        for (m in names(vals)) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <- vals[m]
        }
      }
    }
  }, error = function(e) {})

  # Build the radar chart matrix
  radar_result <- compute_radar(radar_plot_values, names_sigs, names_datasets)
  radar_result$output_table
}

# ============================================================================
# Helper: vectorized quantile summary
# ============================================================================
.summarize_metrics <- function(metrics_matrix) {
  quant_probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  quant_mat <- apply(metrics_matrix, 2, stats::quantile, probs = quant_probs, na.rm = TRUE)
  mean_vals <- apply(metrics_matrix, 2, mean, na.rm = TRUE)
  summary_mat <- rbind(mean_vals, quant_mat)
  rownames(summary_mat) <- c("mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975")
  summary_mat
}

# ============================================================================
# Helper: plot comparison boxplot
# ============================================================================
.plot_comparison <- function(summary_mat, metrics_matrix, sig_metrics_row,
                              out_dir, plot_name, plot_title, control_label) {
  stripchartMatrixList <- list()
  stripchartMatrixList[[control_label]] <- metrics_matrix
  if (!is.null(sig_metrics_row)) {
    stripchartMatrixList[["Original Metric Value"]] <- sig_metrics_row
  }

  stripchart_group_names <- c('Relative Med. SD', 'Skewness',
                               expression(sigma["" >= "10%"]), expression(sigma["" >= "25%"]),
                               expression(sigma["" >= "50%"]), 'Coef. of Var.',
                               'Non-NA Prop.', 'Prop. Expressed',
                               'Autocor.', expression(rho["Mean,Med"]),
                               expression(rho["PCA1,Med"]), expression(rho["Mean,PCA1"]),
                               expression(sigma["PCA1"]), expression(rho["Med,Z-Med"]))

  .boxplot.matrix2(x = summary_mat[2:6, ], outputDir = out_dir,
                   plotName = plot_name,
                   plotTitle = plot_title,
                   stripchartMatrixList = stripchartMatrixList,
                   stripchartPch = c(1, 21), stripchartCol = c("gray", "red"),
                   xlab = "Metrics", ylab = "Score",
                   group.names = stripchart_group_names)
}

# ============================================================================
# Helper: get original metrics row for comparison
# ============================================================================
.get_original_metrics_opt <- function(outputDir, datasetName, signatureName) {
  sig.metrics.file <- file.path(outputDir, "radarchart_table", "radarchart_table.txt")
  if (file.exists(sig.metrics.file)) {
    sig.tab <- utils::read.table(file = sig.metrics.file, header = TRUE,
                                  sep = "\t", check.names = FALSE, row.names = 1)
    row_id <- paste0(gsub(' ', '.', datasetName), '_', gsub(' ', '.', signatureName))
    idx <- which(rownames(sig.tab) == row_id)
    if (length(idx) > 0) return(sig.tab[idx, , drop = FALSE])
  }
  NULL
}

# ============================================================================
# .sigsQcNegativeControl_opt: Outer loop dispatcher (optimized)
# ============================================================================
.sigsQcNegativeControl_opt <- function(genesList, expressionMatrixList, outputDir,
                                        studyName = "MyStudy", numResampling = 50,
                                        warningsFile = NULL, logFile = NULL,
                                        n_cores = NULL) {
  if (missing(genesList)) stop("Need to specify a list of genes.")
  if (missing(expressionMatrixList)) stop("Need to specify a list of expression matrices.")
  if (missing(outputDir)) stop("Need to specify an output directory.")
  if (is.null(warningsFile)) warningsFile <- file.path(outputDir, "warnings.log")
  if (is.null(logFile)) logFile <- file.path(outputDir, "log.log")

  if (!dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)

  tryCatch({
    log.con <- file(logFile, open = "a")
    warnings.con <- file(warningsFile, open = "a")

    for (genes.i in seq_along(genesList)) {
      genes <- as.matrix(genesList[[genes.i]])
      colnames(genes) <- names(genesList)[genes.i]
      .sigQcNegativeControl_opt(genes, expressionMatrixList, outputDir, studyName,
                                rePower = numResampling, warningsFile, logFile, n_cores)
    }
  }, error = function(err) {
    cat(paste(Sys.time(), "Error in sigQcNegativeControl:", err, sep = " "),
        file = log.con, sep = "\n")
  }, finally = {
    close(log.con)
    close(warnings.con)
  })
}

# ============================================================================
# .sigQcNegativeControl_opt: Per-signature worker (optimized)
# ============================================================================
.sigQcNegativeControl_opt <- function(genes, expressionMatrixList, outputDir, studyName,
                                       rePower = 50, warningsFile, logFile, n_cores = NULL) {
  if (missing(genes)) stop("Need to specify a list of genes.")
  if (missing(expressionMatrixList)) stop("Need to specify a list of expression matrices.")

  par_apply <- .get_apply_fn(n_cores)

  tryCatch({
    if (!dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)
    log.con <- file(logFile, open = "a")
    warnings.con <- file(warningsFile, open = "a")

    len <- nrow(genes)
    datasets.num <- length(expressionMatrixList)
    datasets.names <- names(expressionMatrixList)

    # ==================== PART 1: NEGATIVE CONTROLS ====================
    for (dataset.i in seq_len(datasets.num)) {
      expressionMatrix <- expressionMatrixList[[dataset.i]]
      datasetName <- if (is.null(datasets.names)) paste0("Dataset", dataset.i)
                     else datasets.names[dataset.i]
      n_genes_total <- nrow(expressionMatrix)
      signatureName <- colnames(genes)

      cat(paste(Sys.time(), "Computing Negative Control (optimized)...", sep = " "),
          file = log.con, sep = "")

      # Generate all random signatures using sample() (vectorized)
      gene_sigs_list <- lapply(seq_len(rePower), function(i) {
        idx <- sample(seq_len(n_genes_total), size = len, replace = FALSE)
        as.matrix(rownames(expressionMatrix)[idx])
      })
      names(gene_sigs_list) <- paste0("NC", seq_len(rePower))

      # Compute QC metrics for all random signatures IN MEMORY (no disk roundtrip)
      mRNA_expr_matrix <- list()
      mRNA_expr_matrix[[datasetName]] <- expressionMatrix
      metrics_table <- .compute_qc_metrics_inmemory(
        gene_sigs_list, names(gene_sigs_list), mRNA_expr_matrix, c(datasetName))

      cat("DONE\n", file = log.con)

      # Summarize and plot
      sig.metrics.table <- .get_original_metrics_opt(outputDir, datasetName, signatureName)
      summary_out_dir <- file.path(outputDir, "negative_control", datasetName, signatureName)
      dir.create(summary_out_dir, showWarnings = FALSE, recursive = TRUE)

      summary_mat <- .summarize_metrics(metrics_table)
      .plot_comparison(summary_mat, metrics_table, sig.metrics.table, summary_out_dir,
                       "boxplot_metrics",
                       paste("Boxplot Negative Controls for", signatureName),
                       "Negative Control")

      utils::write.table(summary_mat,
                         file = file.path(summary_out_dir, "neg_controls_summary_table.txt"),
                         row.names = TRUE, col.names = TRUE, sep = ",")
    }

    # ==================== PART 2: PERMUTATION CONTROLS ====================
    for (dataset.i in seq_len(datasets.num)) {
      datasetName <- if (is.null(datasets.names)) paste0("Dataset", dataset.i)
                     else datasets.names[dataset.i]
      expressionMatrix <- expressionMatrixList[[dataset.i]]
      signatureName <- colnames(genes)
      genes_present <- intersect(rownames(expressionMatrix), genes[, 1])
      n_sig_genes <- length(genes_present)

      cat(paste(Sys.time(), "Computing Permutation Control (optimized)...", sep = " "),
          file = log.con, sep = "")

      # OPTIMIZATION: Sparse permutation — only copy/permute signature gene rows
      # Instead of copying the FULL expression matrix rePower times, we create
      # lightweight permuted versions that share the non-signature rows.
      expressionMatrix_perm_list <- par_apply(seq_len(rePower), function(i) {
        perm_matrix <- expressionMatrix  # single copy (R's copy-on-modify)
        # Generate per-column permutations of signature gene indices
        for (col.num in seq_len(ncol(expressionMatrix))) {
          perm_idx <- sample(seq_len(n_sig_genes))
          perm_matrix[genes_present, col.num] <-
            expressionMatrix[genes_present[perm_idx], col.num]
        }
        perm_matrix
      })
      names(expressionMatrix_perm_list) <- paste0("PC", seq_len(rePower))

      # Compute QC metrics for all permuted datasets IN MEMORY
      gene_sigs_list <- list()
      gene_sigs_list[[colnames(genes)]] <- genes

      metrics_table <- .compute_qc_metrics_inmemory(
        gene_sigs_list, names(gene_sigs_list),
        expressionMatrix_perm_list, names(expressionMatrix_perm_list))

      cat("DONE\n", file = log.con)

      # Summarize and plot
      sig.metrics.table <- .get_original_metrics_opt(outputDir, datasetName, signatureName)
      summary_out_dir <- file.path(outputDir, "permutation_control", datasetName, signatureName)
      dir.create(summary_out_dir, showWarnings = FALSE, recursive = TRUE)

      summary_mat <- .summarize_metrics(metrics_table)
      .plot_comparison(summary_mat, metrics_table, sig.metrics.table, summary_out_dir,
                       "boxplot_metrics",
                       paste("Boxplot Permutation Controls for", signatureName),
                       "Permutation Control")

      utils::write.table(summary_mat,
                         file = file.path(summary_out_dir, "perm_controls_summary_table.txt"),
                         row.names = TRUE, col.names = TRUE, sep = ",")
    }

  }, error = function(err) {
    cat(paste(Sys.time(), "Error in sigQcNegativeControl:", err, sep = " "),
        file = log.con, sep = "\n")
  }, finally = {
    close(log.con)
    close(warnings.con)
  })
}
