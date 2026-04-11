# test-negative_control.R
# Unit tests for helper functions in R_refactored/sigsQcNegativeControl.R
# Full end-to-end negative control tests are expensive (50 resamplings),
# so we test the building blocks.

fixture <- readRDS(file.path("..", "fixtures", "fixture_small.rds"))

test_that(".compute_qc_metrics produces radarchart_table file", {
  tmp_dir <- tempfile("qc_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE))

  result <- .compute_qc_metrics(
    gene_sigs_list = fixture$gene_sigs_list,
    names_sigs = fixture$names_sigs,
    mRNA_expr_matrix = fixture$mRNA_expr_matrix,
    names_datasets = fixture$names_datasets,
    out_dir = tmp_dir
  )

  # Should have created the radarchart_table file
  table_path <- file.path(tmp_dir, "radarchart_table", "radarchart_table.txt")
  expect_true(file.exists(table_path))

  # Read and validate the table
  tbl <- utils::read.table(table_path, header = TRUE, sep = "\t",
                            check.names = FALSE, row.names = 1)
  expected_rows <- length(fixture$names_sigs) * length(fixture$names_datasets)
  expect_equal(nrow(tbl), expected_rows)
  expect_equal(ncol(tbl), 14)
})

test_that(".compute_qc_metrics returns nested radar values", {
  tmp_dir <- tempfile("qc_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE))

  result <- .compute_qc_metrics(
    gene_sigs_list = fixture$gene_sigs_list,
    names_sigs = fixture$names_sigs,
    mRNA_expr_matrix = fixture$mRNA_expr_matrix,
    names_datasets = fixture$names_datasets,
    out_dir = tmp_dir
  )

  # Result is a nested list of radar values
  expect_type(result, "list")
  for (sig in fixture$names_sigs) {
    expect_true(sig %in% names(result))
    for (ds in fixture$names_datasets) {
      expect_true(ds %in% names(result[[sig]]))
    }
  }
})

test_that(".get_original_metrics returns NULL for non-existent file", {
  result <- .get_original_metrics(tempdir(), "fake_dataset", "fake_sig")
  expect_null(result)
})

test_that(".get_original_metrics retrieves correct row from table", {
  # Create a mock radarchart_table
  tmp_dir <- tempfile("orig_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE))

  table_dir <- file.path(tmp_dir, "radarchart_table")
  dir.create(table_dir, recursive = TRUE)

  mock_table <- matrix(runif(14), nrow = 1, ncol = 14)
  rownames(mock_table) <- "dataset_A_compact_sig"
  colnames(mock_table) <- paste0("metric_", 1:14)
  utils::write.table(mock_table, file = file.path(table_dir, "radarchart_table.txt"),
                      quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

  result <- .get_original_metrics(tmp_dir, "dataset_A", "compact_sig")
  expect_false(is.null(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 14)
})
