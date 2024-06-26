% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reference_projection_utils.R
\name{transfer_col_data}
\alias{transfer_col_data}
\title{Annotate a dataset with the column data from a reference dataset}
\usage{
transfer_col_data(
  ref,
  query,
  columns,
  k = 20,
  ref_reducedDim = "embedding",
  query_reducedDim = "embedding",
  col_data = NULL
)
}
\arguments{
\item{ref}{a matrix or \code{SummarizedExperiment} with existing annotation}

\item{query}{a matrix or \code{SummarizedExperiment} that will be annotated}

\item{columns}{the quoted column names (e.g., \code{vars(class, subclass)}) from
\code{cbind(colData(ref), col_data)} that will be transferred.}

\item{k}{the number of nearest neighbors to consider}

\item{ref_reducedDim, query_reducedDim}{the names of the \code{reducedDim} to pick from
\code{ref} and \code{query}. Only applies if they are \code{SummarizedExperiment}.}

\item{col_data}{additional column annotation for \code{ref}. Default: \code{NULL}}
}
\value{
a tibble with columns specified by the \code{columns} argument
}
\description{
This function finds the \code{k} nearest neighbors from the reference for each observation in the query
and returns the most common row from the \code{columnData(ref)} for each query obseration.
}
