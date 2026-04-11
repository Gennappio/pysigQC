# test-compute_struct.R
# Unit tests for compute_struct() from R_refactored/eval_struct_loc.R
# Note: biclust package is required

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that("compute_struct returns correct structure", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result, "list")
  expect_named(result, c("sig_scores_all_mats", "all_row_names",
                          "biclust_results", "any_biclusters"))
})

test_that("compute_struct all_row_names is union across datasets", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    gene_sig <- fixture$gene_sigs_list[[sig]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    all_rn <- result$all_row_names[[sig]]

    # Every gene in the union should come from at least one dataset
    for (g in all_rn) {
      found <- FALSE
      for (ds in fixture$names_datasets) {
        if (g %in% rownames(fixture$mRNA_expr_matrix[[ds]])) {
          found <- TRUE
          break
        }
      }
      expect_true(found, info = paste("Gene", g, "not found in any dataset for", sig))
    }
  }
})

test_that("compute_struct padded matrices have consistent row count", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    expected_nrow <- length(result$all_row_names[[sig]])
    for (ds in fixture$names_datasets) {
      mat <- result$sig_scores_all_mats[[sig]][[ds]]
      expect_equal(nrow(mat), expected_nrow,
                   info = paste("Row count mismatch for", sig, ds))
    }
  }
})

test_that("compute_struct padded matrices use NA for missing genes (BUG-8 fix)", {
  # Create asymmetric datasets: dataset_X has gene_a,gene_b; dataset_Y has gene_b,gene_c
  # Union for the signature is {gene_a, gene_b, gene_c}
  # dataset_X should have NA for gene_c; dataset_Y should have NA for gene_a
  mat_X <- matrix(rnorm(30), nrow = 3, ncol = 10,
                  dimnames = list(c("gene_a", "gene_b", "gene_d"), paste0("s", 1:10)))
  mat_Y <- matrix(rnorm(30), nrow = 3, ncol = 10,
                  dimnames = list(c("gene_b", "gene_c", "gene_d"), paste0("s", 1:10)))

  test_sigs <- list(test_sig = c("gene_a", "gene_b", "gene_c"))
  test_expr <- list(dataset_X = mat_X, dataset_Y = mat_Y)

  result <- compute_struct(test_sigs, "test_sig", test_expr, c("dataset_X", "dataset_Y"))

  # Union should be {gene_a, gene_b, gene_c}
  expect_equal(length(result$all_row_names[["test_sig"]]), 3)

  # dataset_X: gene_c is missing → should be padded with NA
  mat_padded_X <- result$sig_scores_all_mats[["test_sig"]][["dataset_X"]]
  expect_true("gene_c" %in% rownames(mat_padded_X))
  expect_true(all(is.na(mat_padded_X["gene_c", ])))

  # dataset_Y: gene_a is missing → should be padded with NA
  mat_padded_Y <- result$sig_scores_all_mats[["test_sig"]][["dataset_Y"]]
  expect_true("gene_a" %in% rownames(mat_padded_Y))
  expect_true(all(is.na(mat_padded_Y["gene_a", ])))
})

test_that("compute_struct biclust_results are populated", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    for (ds in fixture$names_datasets) {
      bc <- result$biclust_results[[sig]][[ds]]
      expect_true("z_scores" %in% names(bc))
      expect_true("binarized" %in% names(bc))
      expect_true("biclust_result" %in% names(bc))
      expect_true("threshold" %in% names(bc))
    }
  }
})

test_that("compute_struct biclust z_scores have correct dimensions", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  for (sig in fixture$names_sigs) {
    gene_sig <- fixture$gene_sigs_list[[sig]]
    if (is.matrix(gene_sig)) gene_sig <- as.vector(gene_sig)
    for (ds in fixture$names_datasets) {
      bc <- result$biclust_results[[sig]][[ds]]
      inter <- intersect(gene_sig, rownames(fixture$mRNA_expr_matrix[[ds]]))
      expect_equal(nrow(bc$z_scores), length(inter))
      expect_equal(ncol(bc$z_scores), 10)  # n_samples
    }
  }
})

test_that("compute_struct any_biclusters is logical", {
  result <- compute_struct(fixture$gene_sigs_list, fixture$names_sigs,
                           fixture$mRNA_expr_matrix, fixture$names_datasets)

  expect_type(result$any_biclusters, "logical")
})
