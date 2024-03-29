% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{binVector}
\alias{binVector}
\title{Calculate vector of bins}
\usage{
binVector(x)
}
\arguments{
\item{x}{[\code{arma::vec}]\cr
Vector that should be discretized.}
}
\value{
\code{arma::vec} Vector of discretized x.
}
\description{
This function returns a vector of equally spaced points of length the square root of the size of the vector.
}
\examples{
x = runif(100)
binVector(x)
}
