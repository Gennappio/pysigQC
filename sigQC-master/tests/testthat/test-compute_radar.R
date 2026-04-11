# test-compute_radar.R
# Unit tests for compute_radar() from R_refactored/make_radar_chart_loc.R

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

# Build a complete set of radar_plot_values by running all compute functions
.build_radar_values <- function(fixture) {
  radar_values <- list()
  for (k in seq_along(fixture$names_sigs)) {
    radar_values[[fixture$names_sigs[k]]] <- list()
  }

  var_r <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                       fixture$mRNA_expr_matrix, fixture$names_datasets)
  expr_r <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)
  compact_r <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                   fixture$mRNA_expr_matrix, fixture$names_datasets)
  stan_r <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  # Merge all radar values
  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      all_vals <- c(
        var_r$radar_values[[sig]][[ds]],
        expr_r$radar_values[[sig]][[ds]],
        compact_r$radar_values[[sig]][[ds]],
        stan_r$radar_values[[sig]][[ds]]
      )
      radar_values[[sig]][[ds]] <- all_vals
    }
  }
  radar_values
}

test_that("compute_radar returns correct structure", {
  radar_values <- .build_radar_values(fixture)
  result <- compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)

  expect_type(result, "list")
  expect_true("radar_plot_mat" %in% names(result))
  expect_true("output_table" %in% names(result))
  expect_true("areas" %in% names(result))
  expect_true("legend_labels" %in% names(result))
})

test_that("compute_radar output_table has correct dimensions", {
  radar_values <- .build_radar_values(fixture)
  result <- compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)

  # Should have n_sigs * n_datasets rows, 14 metric columns
  expected_rows <- length(fixture$names_sigs) * length(fixture$names_datasets)
  expect_equal(nrow(result$output_table), expected_rows)
  expect_equal(ncol(result$output_table), 14)
})

test_that("compute_radar values are non-negative (abs applied)", {
  radar_values <- .build_radar_values(fixture)
  result <- compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)

  expect_true(all(result$output_table >= 0, na.rm = TRUE))
})

test_that("compute_radar areas are positive", {
  radar_values <- .build_radar_values(fixture)
  result <- compute_radar(radar_values, fixture$names_sigs, fixture$names_datasets)

  expect_true(all(result$areas >= 0))
})

test_that("compute_radar fills missing metrics with 0", {
  # Give incomplete radar values (missing some metrics)
  radar_values <- list()
  radar_values[["compact_sig"]] <- list()
  radar_values[["compact_sig"]][["dataset_A"]] <- c(sd_median_ratio = 0.5)

  radar_values[["diffuse_sig"]] <- list()
  radar_values[["diffuse_sig"]][["dataset_A"]] <- c(sd_median_ratio = 0.3)

  result <- compute_radar(radar_values, c("compact_sig", "diffuse_sig"), c("dataset_A"))

  # Should have 14 columns (all metrics filled)
  expect_equal(ncol(result$output_table), 14)
})
