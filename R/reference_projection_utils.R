
#' @importFrom glmGamPoi vars
#'
#' @returns see [glmGamPoi::vars].
#'
#' @examples
#'   # `vars` quotes expressions (just like in dplyr)
#'   vars(condition, sample)
#'
#' @export
glmGamPoi::vars


#' Annotate a dataset with the column data from a reference dataset
#'
#' This function finds the `k` nearest neighbors from the reference for each observation in the query
#' and returns the most common row from the `columnData(ref)` for each query obseration.
#'
#' @param ref a matrix or `SummarizedExperiment` with existing annotation
#' @param query a matrix or `SummarizedExperiment` that will be annotated
#' @param columns the quoted column names (e.g., `vars(class, subclass)`) from
#'  `cbind(colData(ref), col_data)` that will be transferred.
#' @param k the number of nearest neighbors to consider
#' @param ref_reducedDim,query_reducedDim the names of the `reducedDim` to pick from
#'   `ref` and `query`. Only applies if they are `SummarizedExperiment`.
#' @param col_data additional column annotation for `ref`. Default: `NULL`
#'
#' @returns a tibble with columns specified by the `columns` argument
#'
#' @export
transfer_col_data <- function(ref, query, columns, k = 20,
                           ref_reducedDim = "embedding", query_reducedDim = "embedding",
                           col_data = NULL){
  query_embedding <- if(is(query, "SummarizedExperiment")){
    t(reducedDim(query, query_reducedDim))
  }else if(is.matrix(query)){
    query
  }else{
    stop("'query' must either be a matrix or a `lemur_fit` object", call. = TRUE)
  }

  if(is(ref, "SummarizedExperiment")){
    col_data <- as.data.frame(glmGamPoi:::get_col_data(ref, col_data))
  }

  id_cols <- lapply(columns, rlang::eval_tidy, data = as.data.frame(col_data))
  if(! all(lengths(id_cols) == 1 | lengths(id_cols) == ncol(ref))){
    stop("The argument 'columns' has lengths ", paste0(lengths(id_cols), collapse = ","), ", which does not match the number of columns ",
         "in 'ref' (", ncol(ref), ")")
  }else if(any(lengths(id_cols) == 1)){
    id_cols[lengths(id_cols) == 1] <- lapply(id_cols[lengths(id_cols) == 1], \(x) rep_len(x, ncol(ref)))
  }
  names(id_cols) <- names(columns) %str_or% purrr::map_chr(columns, rlang::as_label)
  id_cols <- tibble::as_tibble(id_cols)

  ref_embedding <- if(is(ref, "SummarizedExperiment")){
    t(reducedDim(ref, ref_reducedDim))
  }else if(is.matrix(query)){
    ref
  }else{
    stop("'query' must either be a matrix or a `lemur_fit` object")
  }
  stopifnot(nrow(ref_embedding) == nrow(query_embedding))

  index <- BiocNeighbors::buildAnnoy(ref_embedding, transposed = TRUE)
  lookup <- BiocNeighbors::queryAnnoy(precomputed = index, query = query_embedding, transposed = TRUE,
                                      k = k, get.index = TRUE)$index

  row_ids <- vctrs::vec_group_id(id_cols)
  keys <- vctrs::vec_group_loc(id_cols)$key
  transferred_ids <- purrr::map_int(seq_len(nrow(lookup)), \(idx){
    most_common_element(row_ids[lookup[idx,]])
  })
  keys[transferred_ids,]
}


`%str_or%` <- function(x, y){
  if(is.null(x)){
    y
  }else if(is.null(y)){
    x
  }else{
    dplyr::if_else(is.na(x) | x == "", y, x)
  }
}

most_common_element <- function(x){
  ux <- vctrs::vec_unique(x)
  matches <- vctrs::vec_match(x, ux, na_equal = TRUE)
  ux[which.max(table(matches))]
}

