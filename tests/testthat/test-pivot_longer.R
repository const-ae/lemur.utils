test_that("sce_pivot_longer works", {
  n_cells <- 100
  n_genes <- 500
  mat <- matrix(rpois(n = n_cells * n_genes, lambda = 0.3), nrow = n_genes, ncol = n_cells)

  col_data <- data.frame(type = sample(c("A", "B", "C"), size = n_cells, replace = TRUE),
                         cond = sample(c("trt", "ctrl"), size = n_cells, replace = TRUE))
  row_data <- data.frame(symbol = paste0("gene_", seq_len(n_genes)),
                         gene_id = paste0("ENSG000", apply(matrix(sample(0:9, n_genes * 5, TRUE), ncol = 5), 1, paste0, collapse = "")))
  colnames(mat) <- paste0("cell_", seq_len(n_cells))
  rownames(mat) <- paste0("gene_", seq_len(n_genes))

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat, logcounts = log1p(mat)),
                                                    colData = col_data, rowData = row_data)
  tidy_df <- sce_pivot_longer(sce)
  expect_equal(colnames(tidy_df), c("type", "cond"))
  expect_equal(nrow(tidy_df), n_cells)

  tidy_df <- sce_pivot_longer(sce, genes = 1)
  expect_equal(colnames(tidy_df), c("type", "cond", "symbol", "gene_id", "counts"))
  expect_equal(tidy_df$counts, mat[1,], ignore_attr = "names")
  expect_equal(tidy_df[,1:2], tibble::as_tibble(SummarizedExperiment::colData(sce)))


  tidy_df <- sce_pivot_longer(sce, genes = c("gene_3", "gene_5"))
  expect_equal(tidy_df$counts, c(mat[c("gene_3", "gene_5"),]), ignore_attr = "names")

  tidy_df <- sce_pivot_longer(sce, genes = 1:3, cells = lemur::vars(type == "A"))
  expect_equal(tidy_df$counts, c(mat[1:3,sce$type == "A"]), ignore_attr = "names")
  expect_equal(tidy_df$type, rep("A", nrow(tidy_df)))

  tidy_df <- sce_pivot_longer(sce, genes = 1, colnames = "cell_id")
  expect_equal(tidy_df$cell_id, colnames(sce))

  colnames(sce) <- NULL
  tidy_df <- sce_pivot_longer(sce, genes = 1, colnames = "cell_id")
  expect_equal(tidy_df$cell_id, as.character(seq_len(n_cells)))

  tidy_df <- sce_pivot_longer(sce, genes = 1, rownames = "gene_name")
  expect_equal(tidy_df$gene_name, rep("gene_1", n_cells))

  suppressMessages({
    tidy_df <- sce_pivot_longer(sce, genes = 1, rownames = "symbol")
  })
  expect_equal(tidy_df[[3]], tidy_df[[4]])

})
