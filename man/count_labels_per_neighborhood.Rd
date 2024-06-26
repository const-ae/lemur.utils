% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighborhoods_utils.R
\name{count_labels_per_neighborhood}
\alias{count_labels_per_neighborhood}
\title{Count the occurrences of a cell label per neighborhood}
\usage{
count_labels_per_neighborhood(
  data,
  labels,
  fit = NULL,
  cell_names = NULL,
  neighborhood_column_name = "neighborhood",
  id_column_name = "name",
  add_total = TRUE,
  return = c("tidy", "matrix", "sparse_matrix")
)
}
\arguments{
\item{data}{the output of \code{lemur::find_de_neighborhoods}}

\item{labels}{a vector with a label for each cell. It can also be a quoted
column name from \code{colData(fit)} (e.g., \code{vars(cell_type)}).}

\item{fit}{the \code{lemur_fit} object (optional). The column names of fit are used
as the set of possible values in the neighborhood columns.}

\item{cell_names}{a character vector with the set of possible values
in the neighborhood columns.}

\item{neighborhood_column_name, id_column_name}{the identifiers used to extract
the relevant neighborhood and gene name column from \code{data}.}

\item{add_total}{flag indicating if the total counts per label are returned.
If \code{TRUE} and \verb{return="tidy}, an additional column is added; otherwise,
an attribute called "total_counts" is added to matrix output.}

\item{return}{specification what format the data should be returned in.}
}
\description{
Count the occurrences of a cell label per neighborhood
}
