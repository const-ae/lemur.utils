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
  expect_true(tibble::is_tibble(df))
  expect_equal(nrow(df), 10 * 100)
  expect_equal(colnames(df), c("name", "cell", "inside"))

  nei_count <- summarize(df, n = sum(inside), .by = name)$n
  expect_equal(nei_count, lengths(nei$neighborhood))
  inside_cells <- df |>
    mutate(name = as.factor(name)) |>
    dplyr::filter(inside) |>
    summarize(nei = list(cell), .by = name) |>
    tidyr::complete(name, fill = list(nei = list(character(0L))))
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
