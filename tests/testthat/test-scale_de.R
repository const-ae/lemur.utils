test_that("scale_de works", {
  suppressMessages({
    library(ggplot2)
  })
  dat <- tibble::tibble(x = 1:10, y = seq(-3, 5, length.out = 10))
  pl1 <- ggplot(dat, aes(x = x, y = y)) +
    geom_col(aes(fill = y)) +
    scale_fill_de()

  vdiffr::expect_doppelganger("Simple scale_fill_de plot", pl1)

  pl2 <- ggplot(dat, aes(x = x, y = y)) +
    geom_col(aes(fill = y)) +
    scale_colour_de(qlimits = c(0, 0.7))
  vdiffr::expect_doppelganger("scale_colour_de with qlimits", pl2)

  dat <- tibble::tibble(x = 1:10, y = seq(-1, 3, length.out = 10))
  pl3 <- ggplot(dat, aes(x = x, y = exp(y))) +
    geom_col(aes(fill = y)) +
    scale_fill_de(midpoint = 1, qlimits = c(0.5, 1))
  vdiffr::expect_doppelganger("scale_fill_de with midpoint", pl3)
})
