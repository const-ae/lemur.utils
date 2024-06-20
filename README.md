
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lemur.utils

<!-- badges: start -->
<!-- badges: end -->

Helper functions to manage the output of
[`lemur`](https://www.bioconductor.org/packages/lemur/).

## Installation

You can install the development version of lemur.utils from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("const-ae/lemur.utils")
```

# Disclaimer

This package is in an very early stage of development and the API is not
considered stable.

# Documentation

### Neighborhood helpers

#### `neighborhood_to_long_data()`

Convert the neighborhood column from the output of
`lemur::find_de_neighborhoods` to a tidy tibble.

#### `neighborhoods_to_matrix()`

Convert the neighborhood column from the output of
`lemur::find_de_neighborhoods` to a 0/1 matrix.

### Projection onto a reference dataset

#### `transfer_col_data`

After integrating a query and a reference dataset, transfer the
annotation from the reference data to the query data.
