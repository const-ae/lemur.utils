
neighborhoods_intersect <- function(){

}


#' Convert neighborhood column to other format
#'
#' @param data the output of `lemur::find_de_neighborhoods`
#' @param fit the `lemur_fit` object (optional). The column names of fit are used
#'   as the set of possible values in the neighborhood columns.
#' @param cell_names a character vector with the set of possible values
#'   in the neighborhood columns.
#' @param neighborhood_column_name,id_column_name the identifiers used to extract
#'   the relevant neighborhood and gene name column from `data`.
#' @param return_sparse if `TRUE` a `dgCMatrix` is returned. Otherwise a regular
#'   dense matrix is returned.
#' @param only_keep_inside return the data after filtering out the cell
#'   outside the neighborhood for each gene.
#' @param verbose indicator if additional messages are printed.
#'
#' @returns \describe{
#' \item{`neighborhoods_to_matrix()`}{a (sparse) matrix with `length(data[[neighborhood_column_name]])`
#'  rows and one column for each element of `cell_names` (or the set of all cell labels
#'  occurring in `data[[neighborhood_column_name]]`).}
#' \item{`neighborhood_to_long_data()`}{a `tibble` with three columns: _name_, _cell_,
#'  and _inside_.}
#'
#'
#'
#' @export
neighborhoods_to_matrix <- function(data, fit = NULL, cell_names = NULL,
                                    neighborhood_column_name = "neighborhood",
                                    id_column_name = "name",
                                    return_sparse = TRUE,
                                    verbose = TRUE){
  tmp <- get_neigbhoods_and_names_from_df(data, neighborhood_column_name, id_column_name)
  nei <- tmp$neighborhoods
  gene_names <- tmp$names
  cell_names <- get_cell_names(nei, fit, cell_names)


  levels2idx <- seq_along(cell_names)
  names(levels2idx) <- cell_names
  sp_mat <- Matrix::sparseMatrix(i = rep(seq_along(nei), lengths(nei)),
                                 j = levels2idx[unlist(nei)],
                                 x = rep(1, sum(lengths(nei))),
                                 dims = c(length(nei), length(cell_names)),
                                 dimnames = list(gene_names, cell_names),
                                 repr = "C")
  if(return_sparse){
    sp_mat
  }else{
    as.matrix(sp_mat)
  }
}

#' @rdname neighborhoods_to_matrix
neighborhoods_to_long_data <- function(data, fit = NULL, cell_names = NULL,
                                      neighborhood_column_name = "neighborhood",
                                      id_column_name = "name",
                                      only_keep_inside = FALSE,
                                      verbose = TRUE){

  tmp <- get_neigbhoods_and_names_from_df(data, neighborhood_column_name, id_column_name)
  nei <- tmp$neighborhoods
  gene_names <- tmp$names

  if(only_keep_inside){
    res <- list(
      rep(gene_names, lengths(nei)),
      unlist(nei),
      rep(TRUE, sum(lengths(nei)))
    )
  }else{
    cell_names <- get_cell_names(nei, fit, cell_names)
    levels2idx <- seq_along(cell_names)
    names(levels2idx) <- cell_names

    res <- list(
      rep(gene_names, each = length(cell_names)),
      rep(cell_names, times = length(gene_names)),
      rep(FALSE, length(cell_names) * length(gene_names))
    )
    offset <- 0
    for(n in nei){
      res[[3]][offset + levels2idx[n]] <- TRUE
      offset <- offset + length(cell_names)
    }
  }
  names(res) <- c("name", "cell", "inside")
  tibble::as_tibble(res)
}


get_cell_names <- function(neighborhoods, fit, cell_names){
  # Get set of all cells
  all_cell_labels <- unique(unlist(neighborhoods))
  levels <- if(! is.null(cell_names)){
    cell_names
  }else if(! is.null(fit)){
    colnames(fit)
  }else{
    all_cell_labels
  }

  # Check if all cells from nei are inside the levels
  if(! all(all_cell_labels %in% levels)){
    missing_level <- ! all_cell_labels %in% levels
    examples <- paste0(head(all_cell_labels[missing_level], n = 3), collapse = ", ")
    stop(sum(missing_level) , " (", sprintf("%2.0f", mean(missing_level)*100), "%) cells from the neighborhoods ",
         "are not in the cell_names. Here are the first three: \n", examples, call. = FALSE)
  }

  levels
}

get_neigbhoods_and_names_from_df <- function(data, neighborhood_column_name, id_column_name){
  # Get data as a list
  if(! is.data.frame(data) ||
     ! all(c(id_column_name, neighborhood_column_name) %in% colnames(data))){
    stop("data must be a data.frame containing an columns: '",
         neighborhood_column_name, "' and '", id_column_name, "'", call. = FALSE)
  }
  names <- data[[id_column_name]]
  nei <- data[[neighborhood_column_name]]
  if(! is.list(nei)){
    stop("'data[[", neighborhood_column_name, "]]' must be a list.")
  }
  list(neighborhoods = nei, names = names)
}
