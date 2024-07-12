


#' Convert a `lemur_fit` into a tidy datastructure
#'
#' The function takes a `lemur_fit` (or any other `SingleCellExperiment`) and
#' converts it to a tidy tibble.
#'
#' @param sce a `SingleCellExperiment`
#' @param fit the `lemur_fit` produced by `lemur()`
#' @param genes a selection of rows (i.e., genes) from `fit`. It can either
#'   be an index vector, a character vector with `rownames` or a
#'   `dplyr::filter`-like expression wrapped in `vars()` that can refer
#'   to columns of the `rowData(fit)`. Default: `NULL` which means no genes
#'   are included.
#' @param cells a selection of columns (i.e., cells) from `fit`. It can either
#'   be an index vector, a character vector with `colnames` or a
#'   `dplyr::filter`-like expression wrapped in `vars` that can refer to
#'   columns of the `colData(fit)`. Default: `NULL` which means that all
#'   cells are included.
#' @param assays a character vector with the names of assays whose
#'   values are included in the output. To see all options run `assayNames(fit)`.
#' @param reduced_dims a character with the names of the `reducedDim(fit)`
#'  which are included as `matrix` columns in the output. To see all options
#'  run `reducedDimNames(fit)`.
#' @param rownames,colnames string that specifies the column name for
#'  for the `rownames(sce)` and `colnames(sce)` respectively. Default:
#'  `NULL` which means no extra columns are included
#'
#' @returns a `tibble` with the content from
#'   \itemize{
#'     \item `colnames(fit)` if `colnames` is not `NULL`,
#'     \item `colData(fit)` subset according to `cells`,
#'     \item `rownames(fit)` if `rownames` is not `NULL`,
#'     \item `rowData(fit)`subset according to `genes`,
#'     \item one numeric column per entry in `assays`.
#'     \item one matrix-column per entry in `reducedDims`.
#'   }
#'
#'
#' @export
sce_pivot_longer <- function(sce,
                 genes = NULL,
                 cells = NULL,
                 assays = 1,
                 reduced_dims = NULL,
                 rownames = NULL,
                 colnames = NULL){
  col_data <- tibble::as_tibble(SummarizedExperiment::colData(sce))
  col_selection <- if(rlang::is_quosures(cells)){
    evaluate_quosures_to_lgl(cells, col_data, size = ncol(sce))
  }else if(is.null(cells)){
    rep(TRUE, ncol(sce))
  }else{
    cells
  }
  col_selection <- lemur:::convert_subset_to_index(col_selection, colnames(sce))

  skip_genes <- is.null(genes)
  row_data <- tibble::as_tibble(SummarizedExperiment::rowData(sce))
  row_selection <- if(rlang::is_quosures(genes)){
    evaluate_quosures_to_lgl(genes, row_data, size = nrow(sce))
  }else if(is.null(genes)){
    rep(FALSE, nrow(sce))
  }else{
    genes
  }
  row_selection <- lemur:::convert_subset_to_index(row_selection, rownames(sce))
  nrows <- length(row_selection)
  ncols <- length(col_selection)

  if(is.null(reduced_dims)){
    reduced_dim_df <- tibble::tibble(.rows = ncols)
  }else{
    reduced_dim_vals <- lapply(reduced_dims, \(name){
      SingleCellExperiment::reducedDim(sce, name)[col_selection,,drop=FALSE]
    })
    names(reduced_dim_vals) <- if(is.character(reduced_dims)){
      reduced_dims
    }else{
      SingleCellExperiment::reducedDimNames(sce)[reduced_dims]
    }
    reduced_dim_df <- tibble::as_tibble(reduced_dim_vals)
  }

  if(is.null(assays)){
    assay_df <- tibble::tibble(.rows = ncols * nrows)
  }else{
    assay_vals <- lapply(assays, \(name){
      as.vector(SummarizedExperiment::assay(sce, name)[row_selection, col_selection])
    })
    names(assay_vals) <- if(is.character(assays)){
      assays
    }else{
      SummarizedExperiment::assayNames(sce)[assays]
    }
    assay_df <- tibble::as_tibble(assay_vals)
  }

  colnames_df <- if(is.null(colnames)){
    tibble::tibble(.rows = ncols)
  }else if(is.null(colnames(sce))){
    tibble::tibble({{colnames}} := as.character(seq_len(ncols)))
  }else{
    tibble::tibble({{colnames}} := colnames(sce))
  }
  rownames_df <- if(is.null(rownames)){
    tibble::tibble(.rows = nrows)
  }else if(is.null(rownames(sce))){
    tibble::tibble({{rownames}} := as.character(seq_len(nrows)))
  }else{
    tibble::tibble({{rownames}} := rownames(sce))
  }

  if(skip_genes){
    dplyr::bind_cols(
      colnames_df[col_selection,],
      col_data[col_selection,],
      reduced_dim_df
    )
  }else{
    dplyr::bind_cols(
      vctrs::vec_rep_each(colnames_df[col_selection,], nrows),
      vctrs::vec_rep_each(col_data[col_selection,], nrows),
      vctrs::vec_rep_each(reduced_dim_df, nrows),
      vctrs::vec_rep(rownames_df[row_selection,], times = ncols),
      vctrs::vec_rep(row_data[row_selection,], times = ncols),
      tibble::as_tibble(assay_df)
    )
  }

}

#' @rdname sce_pivot_longer
#' @export
fit_pivot_longer <- function(fit,
                             genes = NULL,
                             cells = NULL,
                             assays = fit$use_assay,
                             reduced_dims = "embedding",
                             rownames = NULL,
                             colnames = NULL){
  sce_pivot_longer(fit, genes = genes, cells = cells,
                   assays = assays, reduced_dims = reduced_dims,
                   rownames = rownames, colnames = colnames)
}


evaluate_quosures_to_lgl <- function(filters, data, size){
  MatrixGenerics::colAlls(
    lemur:::mply_dbl(filters, \(f){
      vec <- rlang::eval_tidy(f, data = data)
      if(length(vec) == 1){
        vec <- rep(vec, size)
      }else if(length(vec) != size){
        stop("The conditions inside 'vars(...)' must evaluate to a vector of length 1 or ", size)
      }
      if(! is.logical(vec)){
        stop("The conditions inside 'vars(...)' must evaluate to a logical vector.")
      }
      vec[is.na(vec)] <- FALSE
      vec
    }, ncol = size)
  )
}
