% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rowWiseTensorSparse}
\alias{rowWiseTensorSparse}
\title{Calculating row wise tensor product for sparse matrices}
\usage{
rowWiseTensorSparse(A, B)
}
\arguments{
\item{A}{[\code{matrix}]\cr
First sparse matrix for the tensor product.}

\item{B}{[\code{matrix}]\cr
Second sparse matrix for the tensor product.}
}
\value{
\code{arma::mat} Sparse tensor product
}
\description{
This function calculates the row wise tensor product which can be used
for multivariate smoothing but also for modelling multivariate interactions.
}
\examples{
pen = penaltyMat(10, 2)
}
