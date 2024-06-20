library(dplyr)
library(tidyr)
library(purrr)
library(tibble)


test_that("label counting per neighborhood works", {

  nei <- tibble(name = paste0("gene_", 1:3), neighborhood = list(c(1,2,3,3), c(2), c(4,4,4,4)))
  labels <- rep(letters[1:2], times = 5)
  tidy_res1 <- aggregate_labels_per_neighborhood(nei, labels, return = "tidy")
  mat_res <- aggregate_labels_per_neighborhood(nei, labels, return = "matrix")
  expect_equal(tidy_res1, tibble(name = c("gene_1", "gene_1", "gene_2", "gene_2", "gene_3", "gene_3"),
                                label = c("a", "b", "a", "b", "a", "b"),
                                counts = c(3, 1, 0, 1, 0, 4)))
  expect_equal(colnames(mat_res), c("a", "b"))
  expect_equal(mat_res, matrix(c(3,0,0,1,1,4), ncol = 2), ignore_attr = "dimnames")

  labels <- tibble(a = rep(letters[1:2], times = 5), b = "hello")
  tidy_res2 <- aggregate_labels_per_neighborhood(nei, labels, return = "tidy")
  testthat::expect_warning(
    aggregate_labels_per_neighborhood(nei, labels, return = "matrix")
  )

  nei$neighborhood <- map(nei$neighborhood, \(x) paste0("cell_", x))
  labels <- rep(letters[1:2], times = 5)
  tidy_res3 <- aggregate_labels_per_neighborhood(nei, labels, cell_names = paste0("cell_", 1:10), return = "tidy")
  expect_equal(tidy_res3, tidy_res1)

  nei <- tibble(name ="gene_1", neighborhood = list(c(10)))
  labels <- c('a','b','a')
  expect_error(
    aggregate_labels_per_neighborhood(nei, labels, return = "tidy")
  )

  nei <- make_dummy_neighborhoods_data.frame(nrow = 4, n_cells = 30)
  sce <- SummarizedExperiment::SummarizedExperiment(colData = S4Vectors::DataFrame(a = rep(letters[1:2], length =30)))
  colnames(sce) <- paste0("cell_", 1:30)

  aggregate_labels_per_neighborhood(nei, vars(a, test = "123"), fit = sce)
  mat_res2 <- aggregate_labels_per_neighborhood(nei, vars(a), fit = sce, return = "matrix")
  expect_warning(
    mat_res3 <- aggregate_labels_per_neighborhood(nei, vars(a, test = "123"), fit = sce, return = "matrix")
  )
  expect_equal(mat_res2, mat_res3, ignore_attr = "dimnames")
})

test_that("conversion to matrix works", {
  nei <- make_dummy_neighborhoods_data.frame()
  mat <- neighborhoods_to_matrix(nei, verbose = FALSE)
  expect_true(is(mat, "dgCMatrix"))
  expect_equal(rownames(mat), nei$name)

  mat <- neighborhoods_to_matrix(nei, cell_names = paste0("cell_", 1:100), verbose = FALSE)
  expect_equal(colnames(mat), paste0("cell_", 1:100))

  expect_error(neighborhoods_to_matrix(nei, cell_names = paste0("nonsense_name_", 1:10), verbose = FALSE))

  nei2 <- make_dummy_neighborhoods_data.frame(cell_names = FALSE)
  mat <- neighborhoods_to_matrix(nei2, verbose = FALSE)
  expect_true(all(stringr::str_detect(colnames(mat), "\\d{1,2}")))
})


test_that("conversion to long data", {
  nei <- make_dummy_neighborhoods_data.frame()
  df <- neighborhoods_to_long_data(nei, only_keep_inside = FALSE, verbose = FALSE)
  expect_true(is_tibble(df))
  expect_equal(nrow(df), 10 * 100)
  expect_equal(colnames(df), c("name", "cell", "inside"))

  nei_count <- summarize(df, n = sum(inside), .by = name)$n
  expect_equal(nei_count, lengths(nei$neighborhood))
  inside_cells <- df |>
    mutate(name = as.factor(name)) |>
    filter(inside) |>
    summarize(nei = list(cell), .by = name) |>
    complete(name, fill = list(nei = list(character(0L))))
  for(na in unique(nei$name)){
    expect_setequal(filter(inside_cells, name == na)$nei[[1]], filter(nei, name == na)$neighborhood[[1]])
  }

  manual_mat <- df |>
    mutate(name = factor(name, paste0("gene_", 1:10)),
           cell = factor(cell, paste0("cell_", 1:100))) |>
    arrange(name, cell) |>
    mutate(inside = inside * 1.0) |>
    pivot_wider(id_cols = name, names_from = cell, values_from = "inside") |>
    column_to_rownames("name") |>
    as.matrix()
  expect_equal(manual_mat, neighborhoods_to_matrix(nei, cell_names = paste0("cell_", 1:100),
                                                   verbose = FALSE, return_sparse = FALSE))

  df_filtered <- neighborhoods_to_long_data(nei, only_keep_inside = TRUE, verbose = FALSE)
  df_filtered_manual <- df |> filter(inside)
  expect_equal(arrange(df_filtered, name, cell), arrange(df_filtered_manual, name, cell))
})
