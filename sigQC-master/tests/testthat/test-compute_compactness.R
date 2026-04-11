# test-compute_compactness.R
# Unit tests for compute_compactness() from R_refactored/eval_compactness_loc.R

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("compute_compactness returns correct structure", {
  result <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_named(result, c("radar_values", "autocor_matrices", "rank_product_tables"))
})

test_that("compute_compactness returns autocor_median metric", {
  result <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      expect_true("autocor_median" %in% names(rv))
    }
  }
})

test_that("compact signature has higher autocorrelation than diffuse", {
  result <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (ds in fixture$names_datasets) {
    compact_autocor <- result$radar_values[["compact_sig"]][[ds]]["autocor_median"]
    diffuse_autocor <- result$radar_values[["diffuse_sig"]][[ds]]["autocor_median"]
    # Compact signature should have higher median autocorrelation
    expect_true(compact_autocor > diffuse_autocor,
                info = paste("Dataset:", ds, "compact:", compact_autocor,
                             "diffuse:", diffuse_autocor))
  }
})

test_that("autocorrelation matrices are square and symmetric", {
  result <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      acm <- result$autocor_matrices[[sig]][[ds]]
      expect_equal(nrow(acm), ncol(acm))
      # Diagonal should be 1.0 (self-correlation)
      expect_equal(unname(diag(acm)), rep(1.0, nrow(acm)), tolerance = 1e-10)
      # Symmetric
      expect_equal(acm, t(acm), tolerance = 1e-10)
    }
  }
})

test_that("autocorrelation values are in [-1, 1]", {
  result <- compute_compactness(fixture$gene_sigs_list, fixture$names_sigs,
                                fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      acm <- result$autocor_matrices[[sig]][[ds]]
      expect_true(all(acm >= -1 & acm <= 1, na.rm = TRUE))
    }
  }
})
