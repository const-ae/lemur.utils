

#' Count the occurrences of a cell label per neighborhood
#'
#' @inheritParams neighborhoods_to_matrix
#' @param labels a vector with a label for each cell. It can also be a quoted
#'   column name from `colData(fit)` (e.g., `vars(cell_type)`).
#' @param add_total flag indicating if the total counts per label are returned.
#'   If `TRUE` and `return="tidy`, an additional column is added; otherwise,
#'   an attribute called "total_counts" is added to matrix output.
#' @param return specification what format the data should be returned in.
#'
#'
#' @export
count_labels_per_neighborhood <- function(data, labels, fit = NULL, cell_names = NULL,
                                          neighborhood_column_name = "neighborhood",
                                          id_column_name = "name",
                                          add_total = TRUE,
                                          return = c("tidy", "matrix", "sparse_matrix")){
  tmp <- get_neigbhoods_and_names_from_df(data, neighborhood_column_name, id_column_name)
  neighborhoods <- tmp$neighborhoods
  gene_names <- tmp$names
  cell_names <- get_cell_names(neighborhoods, fit, cell_names)
  stopifnot(vctrs::obj_is_vector(labels))
  return <- rlang::arg_match(return)

  if(! is.list(neighborhoods)){
    stop("'neighborhoods' must be a list where each element is a vector with cell IDs", call. = FALSE)
  }
  if(rlang::is_quosures(labels)){
    col_data <- if(is.null(fit)){
      NULL
    }else{
      as.data.frame(SummarizedExperiment::colData(fit))
    }
    label_names <- purrr::map_chr(labels, rlang::as_label)
    labels <- lapply(labels, rlang::eval_tidy, data = col_data)
    labels <- if(length(labels) == 1){
      labels[[1]]
    }else{
      names(labels) <- label_names
      tibble::as_tibble(labels)
    }
  }

  index_lookup <- is.numeric(vctrs::vec_ptype(neighborhoods[[1]]))
  if(! index_lookup && is.null(cell_names)){
    stop("The neighborhoods contain cell IDs, thus a mapping between cell_id and `labels` is required.", call. = FALSE)
  }else if(! index_lookup && length(cell_names) != vctrs::vec_size(labels)){
    stop("The length of 'labels' and 'cell_names' must match.", call. = FALSE)
  }else if(! index_lookup && length(cell_names) != length(unique(cell_names))){
    stop("The elements in cell_names must be distinct.", call. = FALSE)
  }

  label_ids <- vctrs::vec_group_id(labels)
  n_distinct_ids <- attr(label_ids, "n")
  keys <- vctrs::vec_group_loc(labels)$key
  counts <- if(index_lookup){
    lemur:::mply_dbl(neighborhoods, \(nei){
      if(any(nei <= 0 | nei > vctrs::vec_size(labels))){
        stop("The neighborhood contains an index (", nei[any(nei <= 0 | nei > vctrs::vec_size(labels))][1], ") ",
             "which is outside the bounds of `labels` (1-", vctrs::vec_size(labels), ")", call. = FALSE)
      }
      tabulate(label_ids[nei], n_distinct_ids)
    }, ncol = n_distinct_ids)
  }else{
    lemur:::mply_dbl(neighborhoods, \(nei){
      nei_idx <- vctrs::vec_match(nei, cell_names)
      tabulate(label_ids[nei_idx], n_distinct_ids)
    }, ncol = n_distinct_ids)
  }


  tc <- tabulate(label_ids, n_distinct_ids)
  if(return == "tidy"){
    res <- tibble::tibble(name = vctrs::vec_rep_each(gene_names, ncol(counts)),
           label = vctrs::vec_rep(keys, nrow(counts)),
           counts = c(t(counts)))
    if(add_total){
      res$total_counts <- vctrs::vec_rep(tc, nrow(counts))
    }
    res
  }else{
    if(return == "sparse_matrix"){
      counts <- as(counts, "dgCMatrix")
    }
    if(! is.list(keys) && !is.data.frame(keys)){
      col_names <- purrr::map_chr(seq_len(vctrs::vec_size(keys)), \(idx){
        format(vctrs::vec_slice(keys, idx))
      })
      colnames(counts) <- keys
    }else{
      warning("'keys' is a complex object which cannot be converted to column names. Please use 'return=\"tidy\"")
      colnames(counts) <- paste0("V", seq_len(n_distinct_ids))
    }
    if(add_total){
      names(tc) <- colnames(counts)
      attr(counts, "total_counts") <- tc
    }
    counts
  }
}


#' Convert neighborhood column to other format
#'
#' This function converts the output of [`lemur::find_de_neighborhoods()`]
#' to a format that can be more convenient.
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
#' @param as_factor convert `name` and `cell` columns to factors.
#' @param verbose indicator if additional messages are printed.
#'
#' @returns
#' \describe{
#' \item{`neighborhoods_to_matrix()`}{a (sparse) matrix with `length(data[[neighborhood_column_name]])`
#'  rows and one column for each element of `cell_names` (or the set of all cell labels
#'  occurring in `data[[neighborhood_column_name]]`).}
#' \item{`neighborhood_to_long_data()`}{a `tibble` with three columns: _name_, _cell_,
#'  and _inside_.}
#' }
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
#' @export
neighborhoods_to_long_data <- function(data, fit = NULL, cell_names = NULL,
                                      neighborhood_column_name = "neighborhood",
                                      id_column_name = "name",
                                      only_keep_inside = FALSE,
                                      as_factor = TRUE,
                                      verbose = TRUE){

  tmp <- get_neigbhoods_and_names_from_df(data, neighborhood_column_name, id_column_name)
  nei <- tmp$neighborhoods
  gene_names <- tmp$names
  cell_names <- get_cell_names(nei, fit, cell_names)

  if(only_keep_inside){
    res <- list(
      rep(gene_names, lengths(nei)),
      unlist(nei),
      rep(TRUE, sum(lengths(nei)))
    )
  }else{
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
  if(as_factor){
    res[["name"]] <- factor(res[["name"]], levels = gene_names)
    res[["cell"]] <- factor(res[["cell"]], levels = cell_names)
  }
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
