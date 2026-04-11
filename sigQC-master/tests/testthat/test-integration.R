# test-integration.R
# Integration tests: run the full compute pipeline and verify consistency
# across all modules.

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("full pipeline produces 14 radar metrics per sig-dataset pair", {
  # Run all compute modules
  var_r <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                       fixture$mRNA_expr_matrix, fixture$names_datasets)
  expr_r <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)
  compact_r <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                   fixture$mRNA_expr_matrix, fixture$names_datasets)
  metrics_r <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                               fixture$mRNA_expr_matrix, fixture$names_datasets)
  stan_r <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  # Assemble radar values
  radar_values <- list()
  for (sig in fixture$names_sigs) {
    radar_values[[sig]] <- list()
    for (ds in fixture$names_datasets) {
      all_vals <- c(
        var_r$radar_values[[sig]][[ds]],
        expr_r$radar_values[[sig]][[ds]],
        compact_r$radar_values[[sig]][[ds]],
        metrics_r$radar_values[[sig]][[ds]],
        stan_r$radar_values[[sig]][[ds]]
      )
      radar_values[[sig]][[ds]] <- all_vals
    }
  }

  # Each sig-dataset pair should have 14 metrics:
  # 6 (var) + 2 (expr) + 1 (compact) + 4 (metrics) + 1 (stan) = 14
  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      expect_equal(length(radar_values[[sig]][[ds]]), 14,
                   info = paste("Expected 14 metrics for", sig, ds,
                                "got", length(radar_values[[sig]][[ds]])))
    }
  }
})

test_that("radar chart can be built from full pipeline output", {
  var_r <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                       fixture$mRNA_expr_matrix, fixture$names_datasets)
  expr_r <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)
  compact_r <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                   fixture$mRNA_expr_matrix, fixture$names_datasets)
  metrics_r <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                               fixture$mRNA_expr_matrix, fixture$names_datasets)
  stan_r <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

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

  result <- compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)

  expect_equal(nrow(result$output_table), 4)  # 2 sigs * 2 datasets
  expect_equal(ncol(result$output_table), 14)
  expect_true(all(result$output_table >= 0, na.rm = TRUE))
  expect_length(result$areas, 4)
  expect_true(all(result$areas >= 0))
})

test_that("compute_struct runs on same fixture without error", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_type(result$any_biclusters, "logical")
})

test_that(".compute_qc_metrics produces consistent output with direct pipeline", {
  tmp_dir <- tempfile("integration_")
  on.exit(unlink(tmp_dir, recursive = TRUE))

  # Run via .compute_qc_metrics (as negative control would)
  qc_result <- .compute_qc_metrics(
    gene_sigs_list = fixture$gene_sigs_list,
    names_sigs = fixture$names_sigs,
    mRNA_expr_matrix = fixture$mRNA_expr_matrix,
    names_datasets = fixture$names_datasets,
    out_dir = tmp_dir
  )

  # Read the table it wrote
  table_path <- file.path(tmp_dir, "radarchart_table", "radarchart_table.txt")
  tbl <- utils::read.table(table_path, header = TRUE, sep = "\t",
                            check.names = FALSE, row.names = 1)

  # All values should be non-negative (abs applied in compute_radar)
  expect_true(all(tbl >= 0, na.rm = TRUE))
  expect_equal(ncol(tbl), 14)
})

test_that("pipeline is deterministic across runs", {
  run_pipeline <- function() {
    var_r <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)
    expr_r <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)
    compact_r <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                     fixture$mRNA_expr_matrix, fixture$names_datasets)
    stan_r <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

    radar_values <- list()
    for (sig in fixture$names_sigs) {
      radar_values[[sig]] <- list()
      for (ds in fixture$names_datasets) {
        radar_values[[sig]][[ds]] <- c(
          var_r$radar_values[[sig]][[ds]],
          expr_r$radar_values[[sig]][[ds]],
          compact_r$radar_values[[sig]][[ds]],
          stan_r$radar_values[[sig]][[ds]]
        )
      }
    }
    compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)
  }

  r1 <- run_pipeline()
  r2 <- run_pipeline()

  expect_equal(r1$output_table, r2$output_table)
  expect_equal(r1$areas, r2$areas)
})
