# sigsQcNegativeControl.R — OPTIMIZED
#
# Optimized version of negative and permutation controls, layered on R_refactored/.
# Key optimizations over the refactored version:
#   1. No file I/O roundtrip — QC metrics computed in memory (no radarchart_table.txt
#      write/read between compute and summarize)
#   2. Sparse permutation — only permute signature gene rows, not full matrix copies
#   3. Vectorized quantile computation — single apply() instead of 6 separate calls
#   4. parallel::mclapply() over the compute step, chunked by signature (NC) or
#      dataset (PC) so each worker amortizes per-module overhead across many
#      (sig, dataset) pairs (with lapply fallback on Windows / n_cores <= 1)
#   5. BUG-1 fix: sample() instead of runif()
#   6. Deduplicated summary/plotting logic via helper functions
#
# Compute backend: uses compute_var()/compute_expr()/compute_compactness()/
# compute_stan()/compute_metrics() from R_refactored/. Radar matrix assembly uses
# compute_radar() from R_refactored/make_radar_chart_loc.R (in-memory, no disk write).
#
# Source these before using this file:
#   source("R_refactored/eval_var_loc.R")         # compute_var
#   source("R_refactored/eval_expr_loc.R")        # compute_expr
#   source("R_refactored/eval_compactness_loc.R") # compute_compactness
#   source("R_refactored/eval_stan_loc.R")        # compute_stan
#   source("R_refactored/compare_metrics_loc.R")  # compute_metrics
#   source("R_refactored/make_radar_chart_loc.R") # compute_radar
#   source("R_refactored/boxplot.matrix2.R")      # .boxplot.matrix2
#
# Precondition: .assert_backend() at entry of .sigsQcNegativeControl_opt verifies
# all required compute_* / compute_radar / .boxplot.matrix2 symbols are resolvable
# and fails loudly if the caller forgot to source the backend.

# ============================================================================
# Helper: choose parallel backend. Default is serial (n_cores = 1) because
# end-to-end benchmarks on the medium fixture show that forking overhead and
# per-chunk assembly dominate the compute savings at rePower <= 500. The
# compute step itself scales (~2.5x at 8 cores in isolation), so multi-core
# may still help on much larger rePower / larger expression matrices; callers
# who want to try it pass n_cores explicitly.
# ============================================================================
.get_apply_fn <- function(n_cores = NULL) {
  if (is.null(n_cores) || .Platform$OS.type == "windows") {
    return(list(apply_fn = lapply, n_cores = 1L))
  }
  n_cores <- as.integer(n_cores)
  if (n_cores <= 1L) {
    return(list(apply_fn = lapply, n_cores = 1L))
  }
  list(
    apply_fn = function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = n_cores),
    n_cores = n_cores
  )
}

# ============================================================================
# Helper: merge a compute_*()$radar_values nested list into an accumulator `rv`,
# matching the refactored eval_*_loc merge contract (metric-by-metric assignment).
# ============================================================================
.merge_radar_values <- function(rv, module_radar_values, sig_names, ds_names) {
  for (sn in sig_names) {
    for (dn in ds_names) {
      metrics <- module_radar_values[[sn]][[dn]]
      for (metric in names(metrics)) {
        rv[[sn]][[dn]][metric] <- metrics[metric]
      }
    }
  }
  rv
}

# ============================================================================
# Helper: run all 5 compute_*() modules on a (sub_sigs, sub_ds) partition and
# return a nested radar_plot_values list keyed by sig then dataset. Errors from
# any module propagate up with a stack trace (no silent swallowing).
# ============================================================================
.run_compute_chunk <- function(sub_sigs_list, sub_sig_names,
                                sub_expr_list, sub_ds_names) {
  # Match the refactored merge contract: sig level initialized as list(),
  # dataset level left NULL so metric-by-metric assignment produces a named
  # numeric vector (not a named list, which would break compute_radar).
  rv <- list()
  for (sn in sub_sig_names) rv[[sn]] <- list()

  res_var <- compute_var(sub_sigs_list, sub_sig_names, sub_expr_list, sub_ds_names)
  rv <- .merge_radar_values(rv, res_var$radar_values, sub_sig_names, sub_ds_names)

  res_expr <- compute_expr(sub_sigs_list, sub_sig_names, sub_expr_list, sub_ds_names,
                           thresholds = NULL)
  rv <- .merge_radar_values(rv, res_expr$radar_values, sub_sig_names, sub_ds_names)

  res_comp <- compute_compactness(sub_sigs_list, sub_sig_names, sub_expr_list, sub_ds_names,
                                  logged = TRUE, origin = NULL)
  rv <- .merge_radar_values(rv, res_comp$radar_values, sub_sig_names, sub_ds_names)

  res_met <- compute_metrics(sub_sigs_list, sub_sig_names, sub_expr_list, sub_ds_names,
                             radar_only = TRUE)
  rv <- .merge_radar_values(rv, res_met$radar_values, sub_sig_names, sub_ds_names)

  res_stan <- compute_stan(sub_sigs_list, sub_sig_names, sub_expr_list, sub_ds_names)
  rv <- .merge_radar_values(rv, res_stan$radar_values, sub_sig_names, sub_ds_names)

  rv
}

# ============================================================================
# Helper: precondition guard. Verifies that every backend symbol this file calls
# is resolvable in the current search path, and fails loudly with an actionable
# message otherwise. Call at the entry of the top-level dispatcher.
# ============================================================================
.assert_backend <- function() {
  required <- c("compute_var", "compute_expr", "compute_compactness",
                "compute_stan", "compute_metrics", "compute_radar",
                ".boxplot.matrix2")
  caller_env <- parent.frame()
  missing_fns <- required[!vapply(required, function(x)
    exists(x, envir = caller_env, mode = "function", inherits = TRUE),
    logical(1))]
  if (length(missing_fns) > 0L) {
    stop(sprintf(
      paste0("R_optimized/sigsQcNegativeControl.R: missing backend function(s): %s.\n",
             "Source R_refactored/eval_var_loc.R, eval_expr_loc.R, ",
             "eval_compactness_loc.R, eval_stan_loc.R, compare_metrics_loc.R, ",
             "make_radar_chart_loc.R and boxplot.matrix2.R before calling ",
             ".sigsQcNegativeControl_opt()."),
      paste(missing_fns, collapse = ", ")),
      call. = FALSE)
  }
  invisible(TRUE)
}

# ============================================================================
# Helper: compute QC metrics in memory (no disk I/O), parallelized over the
# larger of (signatures, datasets). Work is split into `n_parallel` chunks so
# each forked worker amortizes per-module overhead across multiple pairs.
# Returns the radar chart output_table directly instead of writing/reading
# radarchart_table.txt.
# ============================================================================
.compute_qc_metrics_inmemory <- function(gene_sigs_list, names_sigs,
                                          mRNA_expr_matrix, names_datasets,
                                          par_apply = lapply, n_parallel = 1L) {
  n_sigs <- length(names_sigs)
  n_ds   <- length(names_datasets)

  # Pick chunking axis: the one with more items. Negative controls have
  # n_sigs == rePower, n_ds == 1; permutation controls have the reverse.
  chunk_on_sigs <- (n_sigs >= n_ds)
  n_axis <- if (chunk_on_sigs) n_sigs else n_ds
  n_chunks <- max(1L, min(as.integer(n_parallel), n_axis))
  groups <- ((seq_len(n_axis) - 1L) %% n_chunks) + 1L
  chunk_idx_list <- split(seq_len(n_axis), groups)

  if (chunk_on_sigs) {
    per_chunk <- par_apply(chunk_idx_list, function(sig_idx) {
      sub_names <- names_sigs[sig_idx]
      .run_compute_chunk(gene_sigs_list[sub_names], sub_names,
                         mRNA_expr_matrix, names_datasets)
    })
  } else {
    per_chunk <- par_apply(chunk_idx_list, function(ds_idx) {
      sub_names <- names_datasets[ds_idx]
      .run_compute_chunk(gene_sigs_list, names_sigs,
                         mRNA_expr_matrix[sub_names], sub_names)
    })
  }

  # Assemble nested radar_plot_values from the per-chunk results
  radar_plot_values <- list()
  for (sn in names_sigs) radar_plot_values[[sn]] <- list()
  for (rv in per_chunk) {
    for (sn in names(rv)) {
      for (dn in names(rv[[sn]])) {
        radar_plot_values[[sn]][[dn]] <- rv[[sn]][[dn]]
      }
    }
  }

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
  .assert_backend()

  stopifnot(
    !missing(genesList), is.list(genesList), length(genesList) > 0L,
    !missing(expressionMatrixList), is.list(expressionMatrixList),
    length(expressionMatrixList) > 0L,
    !missing(outputDir), is.character(outputDir), length(outputDir) == 1L,
    is.numeric(numResampling), length(numResampling) == 1L, numResampling >= 1L
  )
  if (is.null(warningsFile)) warningsFile <- file.path(outputDir, "warnings.log")
  if (is.null(logFile)) logFile <- file.path(outputDir, "log.log")

  if (!dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)

  log.con <- file(logFile, open = "a")
  warnings.con <- file(warningsFile, open = "a")
  on.exit({ try(close(log.con), silent = TRUE)
            try(close(warnings.con), silent = TRUE) }, add = TRUE)

  for (genes.i in seq_along(genesList)) {
    genes <- as.matrix(genesList[[genes.i]])
    colnames(genes) <- names(genesList)[genes.i]
    .sigQcNegativeControl_opt(genes, expressionMatrixList, outputDir, studyName,
                              rePower = numResampling, warningsFile, logFile, n_cores)
  }
}

# ============================================================================
# .sigQcNegativeControl_opt: Per-signature worker (optimized)
# ============================================================================
.sigQcNegativeControl_opt <- function(genes, expressionMatrixList, outputDir, studyName,
                                       rePower = 50, warningsFile, logFile, n_cores = NULL) {
  stopifnot(
    !missing(genes), is.matrix(genes) || is.data.frame(genes), nrow(genes) >= 2L,
    !missing(expressionMatrixList), is.list(expressionMatrixList),
    length(expressionMatrixList) > 0L,
    is.numeric(rePower), length(rePower) == 1L, rePower >= 1L
  )

  apply_spec <- .get_apply_fn(n_cores)
  par_apply  <- apply_spec$apply_fn
  n_parallel <- apply_spec$n_cores

  if (!dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = FALSE, recursive = TRUE)
  log.con <- file(logFile, open = "a")
  warnings.con <- file(warningsFile, open = "a")
  on.exit({ try(close(log.con), silent = TRUE)
            try(close(warnings.con), silent = TRUE) }, add = TRUE)

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

      # Compute QC metrics for all random signatures IN MEMORY (no disk roundtrip).
      # Parallel chunking is applied inside .compute_qc_metrics_inmemory: for the
      # negative-control path this splits signatures across workers.
      mRNA_expr_matrix <- list()
      mRNA_expr_matrix[[datasetName]] <- expressionMatrix
      metrics_table <- .compute_qc_metrics_inmemory(
        gene_sigs_list, names(gene_sigs_list), mRNA_expr_matrix, c(datasetName),
        par_apply = par_apply, n_parallel = n_parallel)

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

      # Sparse permutation — only copy/permute signature gene rows. Run serially
      # because per-iteration work is microseconds: fork overhead of mclapply is
      # a net loss here. The expensive compute step below is the one that runs
      # in parallel.
      expressionMatrix_perm_list <- lapply(seq_len(rePower), function(i) {
        perm_matrix <- expressionMatrix  # single copy (R's copy-on-modify)
        for (col.num in seq_len(ncol(expressionMatrix))) {
          perm_idx <- sample(seq_len(n_sig_genes))
          perm_matrix[genes_present, col.num] <-
            expressionMatrix[genes_present[perm_idx], col.num]
        }
        perm_matrix
      })
      names(expressionMatrix_perm_list) <- paste0("PC", seq_len(rePower))

      # Compute QC metrics for all permuted datasets IN MEMORY. Parallel chunking
      # inside .compute_qc_metrics_inmemory splits datasets across workers here
      # (n_sigs == 1, n_ds == rePower).
      gene_sigs_list <- list()
      gene_sigs_list[[colnames(genes)]] <- genes

      metrics_table <- .compute_qc_metrics_inmemory(
        gene_sigs_list, names(gene_sigs_list),
        expressionMatrix_perm_list, names(expressionMatrix_perm_list),
        par_apply = par_apply, n_parallel = n_parallel)

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
}
