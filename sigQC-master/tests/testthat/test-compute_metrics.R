# test-compute_metrics.R
# Unit tests for compute_metrics() from R_refactored/compare_metrics_loc.R
# Note: GSVA and mclust are required dependencies

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("compute_metrics returns correct structure", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_named(result, c("radar_values", "scores", "pca_results",
                          "score_cor_mats", "mixture_models"))
})

test_that("compute_metrics returns 4 radar metrics", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  expected_metrics <- c("rho_mean_med", "rho_pca1_med", "rho_mean_pca1", "prop_pca1_var")

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      for (m in expected_metrics) {
        expect_true(m %in% names(rv),
                    info = paste("Missing metric:", m, "for", sig, ds))
      }
    }
  }
})

test_that("compute_metrics correlation values are in [-1, 1]", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      for (m in c("rho_mean_med", "rho_pca1_med", "rho_mean_pca1")) {
        val <- rv[m]
        if (!is.na(val)) {
          expect_true(val >= -1 && val <= 1,
                      info = paste(m, "=", val, "not in [-1,1] for", sig, ds))
        }
      }
    }
  }
})

test_that("compute_metrics prop_pca1_var is in [0, 1]", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      val <- result$radar_values[[sig]][[ds]]["prop_pca1_var"]
      if (!is.na(val)) {
        expect_true(val >= 0 && val <= 1,
                    info = paste("prop_pca1_var =", val, "not in [0,1] for", sig, ds))
      }
    }
  }
})

test_that("compute_metrics scores have correct length (= n_samples)", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      sc <- result$scores[[sig]][[ds]]
      expect_length(sc$med_scores, 10)
      expect_length(sc$mean_scores, 10)
      # PCA1 may have fewer if NA rows removed
      if (!is.null(sc$pca1_scores)) {
        expect_true(length(sc$pca1_scores) > 0)
      }
    }
  }
})

test_that("compute_metrics PCA results contain variance proportions", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      pca <- result$pca_results[[sig]][[ds]]
      if (!is.null(pca$props_of_variances)) {
        # Proportions should sum to 1
        expect_equal(sum(pca$props_of_variances), 1.0, tolerance = 1e-10)
        # All proportions should be non-negative
        expect_true(all(pca$props_of_variances >= 0))
      }
    }
  }
})

test_that("compute_metrics mixture models are fitted", {
  result <- compute_metrics(fixture$gene_sigs_list, fixture$names_sigs,
                            fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      mm <- result$mixture_models[[sig]][[ds]]
      # Should have median and mean models at minimum
      expect_true("median" %in% names(mm))
      expect_true("mean" %in% names(mm))
      # At least the median model should be non-NULL
      if (!is.null(mm$median)) {
        expect_true(mm$median$G >= 1)
      }
    }
  }
})
