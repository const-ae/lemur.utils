test_that("label transfer works", {

  ref <- matrix(rnorm(10 * 300), nrow = 10, ncol = 300)
  query <- matrix(rnorm(10 * 30), nrow = 10, ncol = 30)

  label <- sample(letters[1:3], size = 300, replace = TRUE)
  res <- transfer_col_data(ref, query, columns = vars(label, "hello"), k = 5)
  expect_equal(colnames(res), c("label", "\"hello\""))

  col_data <- tibble::tibble(x = label, lst = list(a = 1:4, b = "hello", c = "world")[label])
  res2 <- transfer_col_data(ref, query, columns = vars(x, lst), k = 5, col_data = col_data)
  expect_equal(res$label, res2$x)

  res3 <- transfer_col_data(ref, query, columns = vars(label, greeting = "hello"), k = 5)
  expect_equal(colnames(res3), c("label", "greeting"))
})

test_that("most_common_element works", {
  x <- c(1,3,5,2, 3, 5, 3, 3, 2)
  expect_equal(most_common_element(x), 3)
  expect_equal(most_common_element(letters[x]), "c")
  x <- c(1,NA, NA)
  expect_equal(most_common_element(x), NA_real_)
  x <- c(1,2,1,2,3)
  expect_in(most_common_element(x), c(1,2))

  expect_equal(most_common_element(character(0L)), character(0L))

  x <- list(list(5:8), list(1:3), list(1:3), list(1:3), list(5:8), list(4:10))
  most_common_element(x)
})
