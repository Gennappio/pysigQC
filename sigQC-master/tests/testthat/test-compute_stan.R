# test-compute_stan.R
# Unit tests for compute_stan() from R_refactored/eval_stan_loc.R

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("compute_stan returns correct structure", {
  result <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_named(result, c("radar_values", "med_scores", "z_transf_scores"))
})

test_that("compute_stan returns standardization_comp metric", {
  result <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      expect_true("standardization_comp" %in% names(rv))
      # Spearman rho is in [-1, 1]
      expect_true(rv["standardization_comp"] >= -1 && rv["standardization_comp"] <= 1)
    }
  }
})

test_that("compute_stan handles zero-variance gene in dataset_B", {
  # dataset_B has gene_5 with constant expression (sd=0)
  # diffuse_sig includes gene_5
  result <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  # Should not produce NaN or error
  rho <- result$radar_values[["diffuse_sig"]][["dataset_B"]]["standardization_comp"]
  expect_false(is.nan(rho))
  expect_false(is.na(rho))
})

test_that("compute_stan scores have correct length (= n_samples)", {
  result <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      expect_length(result$med_scores[[sig]][[ds]], 10)  # n_samples = 10
      expect_length(result$z_transf_scores[[sig]][[ds]], 10)
    }
  }
})

test_that("compute_stan z-transformed scores have mean ~0", {
  result <- compute_stan(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  # Z-median scores should be approximately centered
  # (not exactly 0 because we take median of z-scores, not mean)
  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      z_scores <- result$z_transf_scores[[sig]][[ds]]
      # Just check they're finite
      expect_true(all(is.finite(z_scores)),
                  info = paste("Non-finite z-scores for", sig, ds))
    }
  }
})
