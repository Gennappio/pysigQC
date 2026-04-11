# test-compute_var.R
# Unit tests for compute_var() from R_refactored/eval_var_loc.R

# Load fixture
fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("compute_var returns correct structure", {
  result <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                        fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_named(result, c("radar_values", "mean_sd_tables", "all_sd", "all_mean", "inter"))

  # Check nesting: [[sig]][[dataset]]
  for (sig in fixture$names_sigs) {
    expect_true(sig %in% names(result$radar_values))
    for (ds in fixture$names_datasets) {
      expect_true(ds %in% names(result$radar_values[[sig]]))
    }
  }
})

test_that("compute_var returns all 6 radar metrics", {
  result <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                        fixture$mRNA_expr_matrix, fixture$names_datasets)

  expected_metrics <- c("sd_median_ratio", "abs_skewness_ratio",
                        "prop_top_10_percent", "prop_top_25_percent",
                        "prop_top_50_percent", "coeff_of_var_ratio")

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      for (m in expected_metrics) {
        expect_true(m %in% names(rv), info = paste("Missing metric:", m, "for", sig, ds))
      }
    }
  }
})

test_that("compute_var radar values are in [0, 1]", {
  result <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                        fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      rv <- result$radar_values[[sig]][[ds]]
      for (m in names(rv)) {
        val <- rv[m]
        if (!is.na(val)) {
          expect_true(val >= 0 && val <= 1,
                      info = paste(m, "=", val, "not in [0,1] for", sig, ds))
        }
      }
    }
  }
})

test_that("compute_var mean_sd_tables have correct dimensions", {
  result <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                        fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      tbl <- result$mean_sd_tables[[sig]][[ds]]
      expect_equal(ncol(tbl), 2)
      expect_equal(colnames(tbl), c("Mean", "SD"))
      # Number of rows = number of genes in intersection
      inter <- result$inter[[sig]][[ds]]
      expect_equal(nrow(tbl), length(inter))
    }
  }
})

test_that("compute_var handles missing genes correctly", {
  result <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                        fixture$mRNA_expr_matrix, fixture$names_datasets)

  # compact_sig has "gene_missing" which is not in any dataset
  inter <- result$inter[["compact_sig"]][["dataset_A"]]
  expect_false("gene_missing" %in% inter)
  expect_equal(length(inter), 4)  # only gene_1 through gene_4
})

test_that("compute_var is deterministic", {
  result1 <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)
  result2 <- compute_var(fixture$gene_sigs_list, fixture$names_sigs,
                         fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      expect_equal(result1$radar_values[[sig]][[ds]],
                   result2$radar_values[[sig]][[ds]])
    }
  }
})
