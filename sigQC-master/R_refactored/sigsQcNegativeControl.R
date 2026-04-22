# sigsQcNegativeControl.R — REFACTORED
#
# Negative and permutation controls for gene signature QC metrics.
# Uses compute_*() functions directly instead of the duplicated _noplots versions.
# The compute_without_plots.R file is no longer needed.

# ============================================================================
# Helper: compute all radar metrics for a given set of signatures and datasets
# ============================================================================
# This replaces the old .compute_without_plots() function by calling the
# refactored compute_*() functions and assembling the radar chart table.
.compute_qc_metrics <- function(gene_sigs_list, names_sigs, mRNA_expr_matrix,
                                names_datasets, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Initialize radar_plot_values accumulator
  radar_plot_values <- list()
  for (k in seq_along(names_sigs)) {
    radar_plot_values[[names_sigs[k]]] <- list()
  }

  # Run each compute module and merge radar values
  tryCatch({
    var_result <- compute_var(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        for (m in names(var_result$radar_values[[names_sigs[k]]][[names_datasets[i]]])) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <-
            var_result$radar_values[[names_sigs[k]]][[names_datasets[i]]][m]
        }
      }
    }
  }, error = function(err) {})

  tryCatch({
    expr_result <- compute_expr(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        for (m in names(expr_result$radar_values[[names_sigs[k]]][[names_datasets[i]]])) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <-
            expr_result$radar_values[[names_sigs[k]]][[names_datasets[i]]][m]
        }
      }
    }
  }, error = function(err) {})

  tryCatch({
    compact_result <- compute_compactness(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        for (m in names(compact_result$radar_values[[names_sigs[k]]][[names_datasets[i]]])) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <-
            compact_result$radar_values[[names_sigs[k]]][[names_datasets[i]]][m]
        }
      }
    }
  }, error = function(err) {})

  tryCatch({
    # radar_only = TRUE: skip GSVA / mclust / score_cor_mats — those outputs are
    # only consumed by plot_metrics, never by the negative/permutation control
    # summarization. Drops this block from ~1.3s to ~0.09s per call.
    metrics_result <- compute_metrics(gene_sigs_list, names_sigs, mRNA_expr_matrix,
                                      names_datasets, radar_only = TRUE)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        for (m in names(metrics_result$radar_values[[names_sigs[k]]][[names_datasets[i]]])) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <-
            metrics_result$radar_values[[names_sigs[k]]][[names_datasets[i]]][m]
        }
      }
    }
  }, error = function(err) {})

  tryCatch({
    stan_result <- compute_stan(gene_sigs_list, names_sigs, mRNA_expr_matrix, names_datasets)
    for (k in seq_along(names_sigs)) {
      for (i in seq_along(names_datasets)) {
        for (m in names(stan_result$radar_values[[names_sigs[k]]][[names_datasets[i]]])) {
          radar_plot_values[[names_sigs[k]]][[names_datasets[i]]][m] <-
            stan_result$radar_values[[names_sigs[k]]][[names_datasets[i]]][m]
        }
      }
    }
  }, error = function(err) {})

  # Compute radar chart table and write to file
  tryCatch({
    radar_result <- compute_radar(radar_plot_values, names_sigs, names_datasets)
    if (!dir.exists(file.path(out_dir, 'radarchart_table'))) {
      dir.create(file.path(out_dir, 'radarchart_table'))
    }
    utils::write.table(radar_result$output_table,
                       file = file.path(out_dir, 'radarchart_table', 'radarchart_table.txt'),
                       quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)
  }, error = function(err) {})

  radar_plot_values
}

# ============================================================================
# Helper: compute summary quantiles and plot boxplot comparison
# ============================================================================
.summarize_and_plot <- function(metrics.table, sig.metrics.table, out_dir,
                                 plot_name, plot_title, control_label) {
  # Vectorized quantile computation (replaces 6 separate apply calls)
  quant_probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  quant_mat <- apply(metrics.table, 2, stats::quantile, probs = quant_probs, na.rm = TRUE)
  mean_vals <- apply(metrics.table, 2, mean, na.rm = TRUE)
  neg.controls.summary <- rbind(mean_vals, quant_mat)
  rownames(neg.controls.summary) <- c("mean", "Q0.025", "Q0.25", "Q0.5", "Q0.75", "Q0.975")

  # Build stripchartMatrixList for boxplot overlay
  stripchartMatrixList <- list()
  stripchartMatrixList[[control_label]] <- metrics.table
  if (!is.null(sig.metrics.table)) {
    stripchartMatrixList[["Original Metric Value"]] <- sig.metrics.table
  }

  # Metric labels for the 14 radar chart axes
  stripchart_group_names <- c('Relative Med. SD', 'Skewness',
                               expression(sigma["" >= "10%"]), expression(sigma["" >= "25%"]),
                               expression(sigma["" >= "50%"]), 'Coef. of Var.',
                               'Non-NA Prop.', 'Prop. Expressed',
                               'Autocor.', expression(rho["Mean,Med"]),
                               expression(rho["PCA1,Med"]), expression(rho["Mean,PCA1"]),
                               expression(sigma["PCA1"]), expression(rho["Med,Z-Med"]))

  .boxplot.matrix2(x = neg.controls.summary[2:6, ], outputDir = out_dir,
                   plotName = plot_name,
                   plotTitle = plot_title,
                   stripchartMatrixList = stripchartMatrixList,
                   stripchartPch = c(1, 21), stripchartCol = c("gray", "red"),
                   xlab = "Metrics", ylab = "Score",
                   group.names = stripchart_group_names)

  neg.controls.summary
}

# ============================================================================
# .get_original_metrics: retrieve original signature metrics for comparison
# ============================================================================
.get_original_metrics <- function(outputDir, datasetName, signatureName) {
  sig.metrics.file.path <- file.path(outputDir, "radarchart_table", "radarchart_table.txt")
  if (file.exists(sig.metrics.file.path)) {
    sig.metrics.table <- utils::read.table(file = sig.metrics.file.path, header = TRUE,
                                            sep = "\t", check.names = FALSE, row.names = 1)
    index <- which(rownames(sig.metrics.table) ==
                     paste0(gsub(' ', '.', datasetName), '_', gsub(' ', '.', signatureName)))
    if (length(index) > 0) {
      return(sig.metrics.table[index, , drop = FALSE])
    }
  }
  NULL
}

# ============================================================================
# .sigsQcNegativeControl: Outer loop dispatcher
# ============================================================================
.sigsQcNegativeControl <- function(genesList, expressionMatrixList, outputDir, studyName,
                                    numResampling = 50, warningsFile, logFile) {
  if (missing(genesList)) stop("Need to specify a list of genes.")
  if (missing(expressionMatrixList)) stop("Need to specify a list of expression matrices.")
  if (missing(outputDir)) stop("Need to specify an output directory")
  if (missing(studyName)) studyName <- "MyStudy"
  if (missing(warningsFile)) warningsFile <- file.path(outputDir, "warnings.log")
  if (missing(logFile)) logFile <- file.path(outputDir, "log.log")

  tryCatch({
    if (!dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)
    log.con <- file(logFile, open = "a")
    warnings.con <- file(warningsFile, open = "a")

    for (genes.i in seq_along(genesList)) {
      genes <- as.matrix(genesList[[genes.i]])
      colnames(genes) <- names(genesList)[genes.i]
      .sigQcNegativeControl(genes, expressionMatrixList, outputDir, studyName,
                            rePower = numResampling, warningsFile, logFile)
    }
  }, error = function(err) {
    cat("", file = log.con, sep = "\n")
    cat(paste(Sys.time(), "Errors occurred in sigQcNegativeControl:", err, sep = " "),
        file = log.con, sep = "\n")
  }, finally = {
    close(log.con)
    close(warnings.con)
  })
}

# ============================================================================
# .sigQcNegativeControl: Per-signature worker
# ============================================================================
.sigQcNegativeControl <- function(genes, expressionMatrixList, outputDir, studyName,
                                   rePower = 50, warningsFile, logFile) {
  if (missing(genes)) stop("Need to specify a list of genes.")
  if (missing(expressionMatrixList)) stop("Need to specify a list of expression matrices.")
  if (missing(outputDir)) stop("Need to specify an output directory")
  if (missing(studyName)) studyName <- "MyStudy"
  if (missing(warningsFile)) warningsFile <- file.path(outputDir, "warnings.log")
  if (missing(logFile)) logFile <- file.path(outputDir, "log.log")

  tryCatch({
    re.power <- rePower
    if (!dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)
    log.con <- file(logFile, open = "a")
    warnings.con <- file(warningsFile, open = "a")

    len <- dim(genes)[1]
    datasets.num <- length(expressionMatrixList)
    datasets.names <- names(expressionMatrixList)

    # ====================
    # PART 1: NEGATIVE CONTROLS (random gene selection)
    # ====================
    for (dataset.i in seq_len(datasets.num)) {
      expressionMatrix <- expressionMatrixList[[dataset.i]]
      datasetName <- if (is.null(datasets.names)) paste0("Dataset", dataset.i)
                     else datasets.names[dataset.i]

      data.matrix.nrows <- nrow(expressionMatrix)
      signatureName <- colnames(genes)

      # Generate random signatures using sample() (BUG-1 fix)
      gene_sigs_list <- list()
      for (i in seq_len(re.power)) {
        random.index.vector <- sample(seq_len(data.matrix.nrows), size = len, replace = FALSE)
        random.genes <- as.matrix(rownames(expressionMatrix)[random.index.vector])
        gene_sigs_list[[paste0("NC", i)]] <- random.genes
      }
      names_sigs <- names(gene_sigs_list)

      cat(paste(Sys.time(), "Computing the Negative Control...", sep = " "),
          file = log.con, sep = "")

      mRNA_expr_matrix <- list()
      mRNA_expr_matrix[[datasetName]] <- expressionMatrix
      sigQC.out_dir <- file.path(outputDir, "negative_control", datasetName, signatureName, "sigQC")
      dir.create(path = sigQC.out_dir, showWarnings = FALSE, recursive = TRUE)

      # Use compute_*() functions via .compute_qc_metrics instead of .compute_without_plots
      .compute_qc_metrics(gene_sigs_list = gene_sigs_list,
                          names_sigs = names_sigs,
                          mRNA_expr_matrix = mRNA_expr_matrix,
                          names_datasets = c(datasetName),
                          out_dir = sigQC.out_dir)

      cat("DONE", file = log.con, sep = "\n")

      # Read results and summarize
      metrics.file.path <- file.path(sigQC.out_dir, "radarchart_table", "radarchart_table.txt")
      metrics.table <- utils::read.table(file = metrics.file.path, header = TRUE,
                                          sep = "\t", check.names = FALSE, row.names = 1)

      sig.metrics.table <- .get_original_metrics(outputDir, datasetName, signatureName)

      summary_out_dir <- file.path(outputDir, "negative_control", datasetName, signatureName)
      neg.controls.summary <- .summarize_and_plot(
        metrics.table, sig.metrics.table, summary_out_dir,
        "boxplot_metrics",
        paste("Boxplot Negative Controls for", signatureName),
        "Negative Control")

      utils::write.table(neg.controls.summary,
                         file = file.path(summary_out_dir, "neg_controls_summary_table.txt"),
                         row.names = TRUE, col.names = TRUE, sep = ",")
    }

    # ====================
    # PART 2: PERMUTATION CONTROLS (gene label shuffling)
    # ====================
    for (dataset.i in seq_len(datasets.num)) {
      datasetName <- if (is.null(datasets.names)) paste0("Dataset", dataset.i)
                     else datasets.names[dataset.i]
      expressionMatrix <- expressionMatrixList[[dataset.i]]
      signatureName <- colnames(genes)

      # Generate permuted expression matrices
      expressionMatrix_perm_list <- list()
      for (i in seq_len(re.power)) {
        expressionMatrix_perm <- expressionMatrix
        genes_present <- intersect(rownames(expressionMatrix), genes[, 1])
        new_ordering <- replicate(sample(seq_along(genes_present)),
                                  n = ncol(expressionMatrix))
        for (col.num in seq_len(ncol(expressionMatrix))) {
          expressionMatrix_perm[genes_present, col.num] <-
            expressionMatrix[genes_present[new_ordering[, col.num]], col.num]
        }
        expressionMatrix_perm_list[[paste0("PC", i)]] <- expressionMatrix_perm
      }

      sigQC.out_dir <- file.path(outputDir, "permutation_control", datasetName, signatureName, "sigQC")
      dir.create(path = sigQC.out_dir, showWarnings = FALSE, recursive = TRUE)

      gene_sigs_list <- list()
      gene_sigs_list[[colnames(genes)]] <- genes

      .compute_qc_metrics(gene_sigs_list = gene_sigs_list,
                          names_sigs = names(gene_sigs_list),
                          mRNA_expr_matrix = expressionMatrix_perm_list,
                          names_datasets = names(expressionMatrix_perm_list),
                          out_dir = sigQC.out_dir)

      cat("DONE", file = log.con, sep = "\n")

      metrics.file.path <- file.path(sigQC.out_dir, "radarchart_table", "radarchart_table.txt")
      metrics.table <- utils::read.table(file = metrics.file.path, header = TRUE,
                                          sep = "\t", check.names = FALSE, row.names = 1)

      sig.metrics.table <- .get_original_metrics(outputDir, datasetName, signatureName)

      summary_out_dir <- file.path(outputDir, "permutation_control", datasetName, signatureName)
      neg.controls.summary <- .summarize_and_plot(
        metrics.table, sig.metrics.table, summary_out_dir,
        "boxplot_metrics",
        paste("Boxplot Permutation Controls for", signatureName),
        "Permutation Control")

      utils::write.table(neg.controls.summary,
                         file = file.path(summary_out_dir, "perm_controls_summary_table.txt"),
                         row.names = TRUE, col.names = TRUE, sep = ",")
    }

  }, error = function(err) {
    cat("", file = log.con, sep = "\n")
    cat(paste(Sys.time(), "Errors occurred in sigQcNegativeControl:", err, sep = " "),
        file = log.con, sep = "\n")
  }, finally = {
    close(log.con)
    close(warnings.con)
  })
}
