# test-compute_expr.R
# Unit tests for compute_expr() from R_refactored/eval_expr_loc.R

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("compute_expr returns correct structure", {
  result <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_named(result, c("radar_values", "na_proportions", "expr_proportions", "thresholds"))
})

test_that("compute_expr returns 2 radar metrics", {
  result <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      expect_true("med_prop_na" %in% names(rv))
      expect_true("med_prop_above_med" %in% names(rv))
    }
  }
})

test_that("compute_expr radar values are in [0, 1]", {
  result <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      for (m in names(rv)) {
        expect_true(rv[m] >= 0 && rv[m] <= 1,
                    info = paste(m, "=", rv[m], "not in [0,1]"))
      }
    }
  }
})

test_that("compute_expr detects NA values in dataset_A", {
  result <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  # compact_sig in dataset_A has gene_3 with 1 NA out of 10 samples
  na_props <- result$na_proportions[["compact_sig"]][["dataset_A"]]
  # gene_missing should be 1.0 (fully NA)
  expect_equal(unname(na_props["gene_missing"]), 1.0)
})

test_that("compute_expr thresholds are computed per dataset", {
  result <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_length(result$thresholds, 2)
  expect_named(result$thresholds, fixture$names_datasets)
  # Thresholds should be positive (expression data is positive)
  expect_true(all(result$thresholds > 0))
})

test_that("compute_expr respects user-provided thresholds", {
  custom_thresholds <- c(3.0, 5.0)
  names(custom_thresholds) <- fixture$names_datasets

  result <- compute_expr(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets,
                         thresholds = custom_thresholds)

  expect_equal(result$thresholds, custom_thresholds)
})
