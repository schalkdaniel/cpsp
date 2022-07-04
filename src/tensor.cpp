#ifndef TENSOR_CPP_
#define TENSOR_CPP_

#include <RcppArmadillo.h>

//' Calculating row wise tensor product
//'
//' This function calculates the row wise tensor product which can be used
//' for multivariate smoothing but also for modelling multivariate interactions.
//'
//' @param A [\code{matrix}]\cr
//'   First matrix for the tensor product.
//' @param B [\code{matrix}]\cr
//'   Second matrix for the tensor product.
//' @return \code{arma::mat} Tensor product
//' @examples
//' pen = penaltyMat(10, 2)
//' @export
// [[Rcpp::export]]
arma::mat rowWiseTensor (const arma::mat& A, const arma::mat& B)
{
  // Variables
  arma::mat out;
  arma::rowvec vecA = arma::rowvec(A.n_cols, arma::fill::ones);
  arma::rowvec vecB = arma::rowvec(B.n_cols, arma::fill::ones);

  // Multiply both kronecker products element-wise
  out = arma::kron(A,vecB) % arma::kron(vecA, B);

  return out;
}

//' Calculating row wise tensor product for sparse matrices
//'
//' This function calculates the row wise tensor product which can be used
//' for multivariate smoothing but also for modelling multivariate interactions.
//'
//' @param A [\code{matrix}]\cr
//'   First sparse matrix for the tensor product.
//' @param B [\code{matrix}]\cr
//'   Second sparse matrix for the tensor product.
//' @return \code{arma::mat} Sparse tensor product
//' @examples
//' pen = penaltyMat(10, 2)
//' @export
// [[Rcpp::export]]
arma::sp_mat rowWiseTensorSparse (const arma::sp_mat& A, const arma::sp_mat& B)
{
  // Variables
  arma::rowvec vecA = arma::rowvec(A.n_cols, arma::fill::ones);
  arma::rowvec vecB = arma::rowvec(B.n_cols, arma::fill::ones);

  arma::sp_mat vecAsparse = arma::sp_mat(vecA);
  arma::sp_mat vecBsparse = arma::sp_mat(vecB);

  // Multiply both kronecker products element-wise
  arma::sp_mat out = arma::kron(A,vecBsparse) % arma::kron(vecAsparse, B);

  return out;
}

# endif // TENSOR_CPP_
