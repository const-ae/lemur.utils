
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

I will demonstrate the functions using the data by Kang et al. (2018).

``` r
library(lemur)
library(lemur.utils)
library(SingleCellExperiment)
library(tidyverse)
# Prepare the data
sce <- muscData::Kang18_8vs8()
logcounts(sce) <- transformGamPoi::shifted_log_transform(sce)
hvg <- order(-rowVars(logcounts(sce)))
sce <- sce[hvg[1:500],]

fit <- lemur(sce, design = ~ stim, n_embedding = 10, verbose = FALSE)
fit <- align_harmony(fit)
fit <- test_de(fit, contrast = cond(stim = "stim") - cond(stim = "ctrl"))
nei <- find_de_neighborhoods(fit, group_by = vars(ind))

as_tibble(nei)
#> # A tibble: 500 × 13
#>    name   neighborhood n_cells sel_statistic     pval adj_pval f_statistic   df1
#>    <chr>  <I<list>>      <int>         <dbl>    <dbl>    <dbl>       <dbl> <int>
#>  1 FTL    <chr>          18155         -120. 2.55e- 3  3.96e-3        12.4     1
#>  2 FTH1   <chr>          20807         -615. 2.11e- 7  9.00e-7        67.4     1
#>  3 ISG15  <chr>           6532          765. 8.02e-10  8.02e-9       142.      1
#>  4 CXCL10 <chr>           8139          234. 8.28e- 7  3.01e-6        55.3     1
#>  5 TIMP1  <chr>           8369         -299. 5.73e- 5  1.36e-4        27.8     1
#>  6 TMSB4X <chr>           8225         -124. 2.18e- 5  5.73e-5        32.9     1
#>  7 CCL2   <chr>           5611          167. 2.21e- 6  7.28e-6        47.7     1
#>  8 RPL3   <chr>          26589         -541. 4.01e- 3  5.99e-3        11.0     1
#>  9 CCL8   <chr>           3786          208. 3.05e- 8  1.69e-7        88.2     1
#> 10 RPS18… <chr>           8980         -209. 1.46e- 3  2.42e-3        14.2     1
#> # ℹ 490 more rows
#> # ℹ 5 more variables: df2 <dbl>, lfc <dbl>, did_pval <dbl>, did_adj_pval <dbl>,
#> #   did_lfc <dbl>
```

### Neighborhood helpers

#### `neighborhoods_to_long_data()`

Convert the neighborhood column from the output of
`lemur::find_de_neighborhoods` to a tidy tibble.

``` r
# By default the long data contains all gene / cell combinations
neighborhoods_to_long_data(nei, fit = fit)
#> # A tibble: 14,532,500 × 3
#>    name  cell             inside
#>    <fct> <fct>            <lgl> 
#>  1 FTL   AAACATACAATGCC-1 FALSE 
#>  2 FTL   AAACATACATTTCC-1 TRUE  
#>  3 FTL   AAACATACCAGAAA-1 TRUE  
#>  4 FTL   AAACATACCAGCTA-1 TRUE  
#>  5 FTL   AAACATACCATGCA-1 TRUE  
#>  6 FTL   AAACATACCTCGCT-1 FALSE 
#>  7 FTL   AAACATACCTGGTA-1 FALSE 
#>  8 FTL   AAACATACGATGAA-1 FALSE 
#>  9 FTL   AAACATACGCCAAT-1 TRUE  
#> 10 FTL   AAACATACGCTTCC-1 TRUE  
#> # ℹ 14,532,490 more rows
# `only_keep_inside` filters out `inside == FALSE` and produces a smaller tibble
neighborhoods_to_long_data(nei, fit = fit, only_keep_inside = TRUE)
#> # A tibble: 6,332,345 × 3
#>    name  cell             inside
#>    <fct> <fct>            <lgl> 
#>  1 FTL   AAACATACATTTCC-1 TRUE  
#>  2 FTL   AAACATACCAGAAA-1 TRUE  
#>  3 FTL   AAACATACCAGCTA-1 TRUE  
#>  4 FTL   AAACATACCATGCA-1 TRUE  
#>  5 FTL   AAACATACGCCAAT-1 TRUE  
#>  6 FTL   AAACATACGCTTCC-1 TRUE  
#>  7 FTL   AAACATACGGCATT-1 TRUE  
#>  8 FTL   AAACATACGTGTAC-1 TRUE  
#>  9 FTL   AAACATACGTTGTG-1 TRUE  
#> 10 FTL   AAACATACTGCGTA-1 TRUE  
#> # ℹ 6,332,335 more rows
```

#### `neighborhoods_to_matrix()`

Convert the neighborhood column from the output of
`lemur::find_de_neighborhoods` to a 0/1 matrix.

``` r
neighborhoods_to_matrix(nei, fit = fit)
#> 500 x 29065 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 33 column names 'AAACATACAATGCC-1', 'AAACATACATTTCC-1', 'AAACATACCAGAAA-1' ... ]]
#>   [[ suppressing 33 column names 'AAACATACAATGCC-1', 'AAACATACATTTCC-1', 'AAACATACCAGAAA-1' ... ]]
#>                                                                              
#> FTL  . 1 1 1 1 . . . 1 1 1 1 1 1 . 1 1 1 1 1 . 1 . 1 1 . 1 1 . . 1 . 1 ......
#> FTH1 . 1 1 1 . 1 1 1 1 1 1 . . . 1 1 . . 1 1 1 . 1 1 1 . 1 . 1 1 . 1 1 ......
#> 
#>  ..............................
#>  ........suppressing 29032 columns and 497 rows in show(); maybe adjust options(max.print=, width=)
#>  ..............................
#>   [[ suppressing 33 column names 'AAACATACAATGCC-1', 'AAACATACATTTCC-1', 'AAACATACCAGAAA-1' ... ]]
#>                                                                                 
#> HNRNPA3 . 1 1 1 1 1 1 . 1 . 1 1 1 1 . . 1 1 1 1 . 1 . . 1 1 1 1 . 1 1 . 1 ......
```

#### `count_labels_per_neighborhood()`

Count the occurrences of a cell label per neighborhood.

``` r
count_labels_per_neighborhood(nei, labels = vars(cell), fit = fit)
#> # A tibble: 4,500 × 4
#>    name  label             counts total_counts
#>    <chr> <fct>              <dbl>        <int>
#>  1 FTL   CD4 T cells         6670        12033
#>  2 FTL   CD14+ Monocytes     5322         6447
#>  3 FTL   Dendritic cells       73          472
#>  4 FTL   NK cells            1400         2330
#>  5 FTL   CD8 T cells         1007         2634
#>  6 FTL   B cells             1961         2880
#>  7 FTL   Megakaryocytes       182          346
#>  8 FTL   FCGR3A+ Monocytes   1531         1914
#>  9 FTL   <NA>                   9            9
#> 10 FTH1  CD4 T cells         8346        12033
#> # ℹ 4,490 more rows
```

### Projection onto a reference dataset

#### `transfer_col_data`

After integrating a query and a reference dataset, transfer the
annotation from the reference data to the query data. The `ref` and
`query` data must contain a shared embedding, produced by manually
integrating them with Seurat, harmony, or LEMUR.

``` r
# This is a completely simulated example to demonstrate how to call the `transfer_col_data` function
ref_sce <- SingleCellExperiment(list(logcounts = matrix(1, nrow = 400, ncol = 300)),
                                colData = DataFrame(celltype = sample(c("A", "B", "C"), size = 300, replace = TRUE),
                                                    origin = sample(c("foo", "bar"), size = 300, replace = TRUE)),
                                reducedDims = list(embedding = t(matrix(rnorm(10 * 300), nrow = 10, ncol = 300))))

transfer_col_data(ref_sce, fit, columns = vars(celltype, origin))
#> # A tibble: 29,065 × 2
#>    celltype origin
#>    <chr>    <chr> 
#>  1 B        bar   
#>  2 A        foo   
#>  3 A        foo   
#>  4 A        foo   
#>  5 B        bar   
#>  6 A        foo   
#>  7 A        foo   
#>  8 A        bar   
#>  9 A        foo   
#> 10 B        bar   
#> # ℹ 29,055 more rows
```

# Session Info

``` r
sessionInfo()
#> R version 4.3.2 (2023-10-31)
#> Platform: x86_64-apple-darwin20 (64-bit)
#> Running under: macOS Sonoma 14.5
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Europe/Berlin
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] muscData_1.16.0      ExperimentHub_2.10.0 AnnotationHub_3.10.0
#>  [4] BiocFileCache_2.10.1 dbplyr_2.4.0         lubridate_1.9.3     
#>  [7] forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4         
#> [10] purrr_1.0.2         
#>  [ reached getOption("max.print") -- omitted 17 entries ]
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.2.2                 bitops_1.0-7             
#>  [3] transformGamPoi_1.6.0     rlang_1.1.3              
#>  [5] magrittr_2.0.3            compiler_4.3.2           
#>  [7] RSQLite_2.3.5             DelayedMatrixStats_1.24.0
#>  [9] png_0.1-8                 vctrs_0.6.5              
#>  [ reached getOption("max.print") -- omitted 75 entries ]
```
