#ifndef BINNING_
#define BINNING_

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>


//' Calculate vector of bins of specific size
//'
//' This function returns a vector of equally spaced points of length n_bins.
//'
//' @param x [\code{arma::vec}]\cr
//'   Vector that should be discretized.
//' @param n_bins [\code{unsigned int}]\cr
//'   Number of unique points for binning the vector x.
//' @return \code{arma::vec} Vector of discretized x.
//' @examples
//' x = runif(100)
//' binVectorCustom(x, 10)
//' @export
// [[Rcpp::export]]
arma::vec binVectorCustom (const arma::vec& x, const unsigned int n_bins)
{
  // TODO: Check if n_bins is set correctly
  return arma::linspace(arma::min(x), arma::max(x), n_bins);
}

//' Calculate vector of bins
//'
//' This function returns a vector of equally spaced points of length the square root of the size of the vector.
//'
//' @param x [\code{arma::vec}]\cr
//'   Vector that should be discretized.
//' @return \code{arma::vec} Vector of discretized x.
//' @examples
//' x = runif(100)
//' binVector(x)
//' @export
// [[Rcpp::export]]
arma::vec binVector (const arma::vec& x)
{
  const unsigned int n_bins = std::floor(std::sqrt(x.size()));
  return binVectorCustom(x, n_bins);
}

//' Calculate index vector for binned vector
//'
//' This function returns the indexes of the unique values to the complete binned vector.
//'
//' @param x [\code{arma::vec}]\cr
//'   Vector that should be discretized.
//' @param x_bins [\code{arma::vec}]\cr
//'   Vector of unique values for binning.
//' @return \code{arma::uvec} Index vector.
//' @examples
//' x = runif(10000)
//' bins = binVector(x)
//' idx = calculateIndexVector(x, bins)
//' head(data.frame(x = x, bins = bins[idx + 1]))
//' @export
// [[Rcpp::export]]
arma::uvec calculateIndexVector (const arma::vec& x, const arma::vec& x_bins)
{
  arma::uvec idx(x.size(), arma::fill::zeros);
  const double delta = (x_bins(1) - x_bins(0)) / 2;

  for (unsigned int i = 0; i < x.size(); i++) {
    unsigned int j = 0;
    while ((x_bins(j) + delta) < x(i)) { j += 1; }
    idx(i) = j;
  }
  return idx;
}

#endif // BINNED_
