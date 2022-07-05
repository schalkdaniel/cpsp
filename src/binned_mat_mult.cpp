#ifndef BINNED_MAT_MULT_
#define BINNED_MAT_MULT_

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

//' Calculating binned matrix product
//'
//' This function calculates the matrix product using Algorithm 3 of Zheyuan Li, Simon N. Wood: "Faster
//' model matrix crossproducts for large generalized linear models with discretized covariates". The idea
//' is to compute just on the unique rows of X by also using an index vector to map to the original matrix.
//' The algorithm implemented here is a small adaption of the original algorithm. Instead of calculating $XW$
//' which again, needs to be transposed, we directly calculate $X^TW$ to avoid another transposing step.
//'
//' @param X [\code{arma::mat}]\cr
//'   Matrix X.
//' @param k [\code{arma::uvec}]\cr
//'   Index vector for mapping to original matrix $X_o(i,) = X(k(i),.)$.
//' @param w [\code{arma::vec}]\cr
//'   Vector of weights that are accumulated.
//' @param use_fast_acc [\code{bool}]\cr
//'   Flag to indicate whether to use the original or adopted algorithm.
//' @return \code{arma::mat} Matrix Product $X^TWX$.
//' @examples
//' nsim = 1e6L
//' nunique = trunc(sqrt(nsim))
//'
//' xunique = runif(n = nunique, min = 0, max = 10)
//' k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)
//'
//' X = poly(x = xunique, degree = 20L)
//'
//' binnedMatMult(X = X, k = k-1, w = 1)
//' @export
// [[Rcpp::export]]
arma::mat binnedMatMult (const arma::mat& X, const arma::uvec& k, const arma::vec& w, const bool use_fast_acc = false)
{
  unsigned int n = k.size();
  unsigned int ind;

  if (use_fast_acc) {
    arma::colvec wcum(X.n_rows, arma::fill::zeros);
    if ( (w.size() == 1) && (w(0) == 1) ) {
      for (unsigned int i = 0; i < n; i++) {
        ind = k(i);
        wcum(ind) += 1;
      }
    } else {
      for (unsigned int i = 0; i < n; i++) {
        ind = k(i);
        wcum(ind) += w(i);
      }
    }
    return arma::trans(X.each_col() % wcum) * X;
  }


  arma::mat L(X.n_cols, X.n_rows, arma::fill::zeros);
  if ( (w.size() == 1) && (w(0) == 1) ) {
    // std::cout << "Weight of size 1:" << std::endl;
    for (unsigned int i = 0; i < n; i++) {
       ind = k(i);
       L.col(ind) += arma::trans(X.row(ind));
    }
  } else {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      L.col(ind) += arma::trans(X.row(ind)) * w(i);
    }
  }
  return L * X;
}

//' Calculating binned matrix product for response term
//'
//' This function calculates the matrix product using Algorithm 3 of Zheyuan Li, Simon N. Wood: "Faster
//' model matrix crossproducts for large generalized linear models with discretized covariates". The idea
//' is to compute just on the unique rows of X by also using an index vector to map to the original matrix.
//' The algorithm implemented here is a small adaption of the original algorithm. Instead of calculating $XW$
//' which again, needs to be transposed, we directly calculate $X^TW$ to avoid another transposing step. In addition
//' to the original algorithm the algorithm here directly calculates the crossproduct with the response.
//'
//' @param X [\code{arma::mat}]\cr
//'   Matrix X.
//' @param y [\code{arma::vec}]\cr
//'   Response vector y.
//' @param k [\code{arma::uvec}]\cr
//'   Index vector for mapping to original matrix $X_o(i,) = X(k(i),.)$.
//' @param w [\code{arma::vec}]\cr
//'   Vector of weights that are accumulated.
//' @return \code{arma::mat} Matrix Product $X^TWX$.
//' @examples
//' nsim = 1e6L
//' nunique = trunc(sqrt(nsim))
//'
//' xunique = runif(n = nunique, min = 0, max = 10)
//' k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)
//'
//' X = poly(x = xunique, degree = 20L)
//' y = runif(nsim)
//'
//' binnedMatMultResponse(X = X, y = y, k = k-1, w = 1)
//' @export
// [[Rcpp::export]]
arma::mat binnedMatMultResponse (const arma::mat& X, const arma::vec& y,  const arma::uvec& k, const arma::vec& w)
{
  unsigned int n = k.size();
  unsigned int ind;

  arma::rowvec wcum(X.n_rows, arma::fill::zeros);

  if ( (w.size() == 1) && (w(0) == 1) ) {
    // std::cout << "Weight of size 1:" << std::endl;
    for (unsigned int i = 0; i < n; i++) {
       ind = k(i);
       wcum(ind) += y(i);
    }
  } else {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      wcum(ind) += w(i) * y(i);
    }
  }
  return wcum * X;
}




//' Calculating sparse binned matrix product
//'
//' This function calculates the matrix product (for sparse matrizes) using Algorithm 3 of Zheyuan Li, Simon N. Wood: "Faster
//' model matrix crossproducts for large generalized linear models with discretized covariates". The idea
//' is to compute just on the unique rows of X by also using an index vector to map to the original matrix.
//' The algorithm implemented here is a small adaption of the original algorithm. Instead of calculating $XW$
//' which again, needs to be transposed, we directly calculate $X^TW$ to avoid another transposing step.
//'
//' @param X [\code{arma::sp_mat}]\cr
//'   Matrix X.
//' @param k [\code{arma::uvec}]\cr
//'   Index vector for mapping to original matrix $X_o(i,) = X(k(i),.)$.
//' @param w [\code{arma::vec}]\cr
//'   Vector of weights that are accumulated.
//' @return \code{arma::mat} Matrix Product $X^TWX$.
//' @examples
//' nsim = 1e6L
//' nunique = trunc(sqrt(nsim))
//'
//' xunique = runif(n = nunique, min = 0, max = 10)
//' k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)
//'
//' X = poly(x = xunique, degree = 20L)
//'
//' binnedMatMult(X = X, k = k-1, w = 1)
//' @export
// [[Rcpp::export]]
arma::mat binnedSparseMatMult (const arma::sp_mat& X, const arma::uvec& k, const arma::vec& w)
{
  const unsigned int n = k.size();
  const unsigned int n_unique = X.n_cols;
  unsigned int ind;

  arma::sp_mat sp_out(X);
  arma::colvec wcum(n_unique, arma::fill::zeros);

  if ( (w.size() == 1) && (w(0) == 1) ) {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      wcum(ind) += 1;
    }
  } else {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      wcum(ind) += w(i);
    }
  }
  for (unsigned int i = 0; i < n_unique; i++) {
    sp_out.col(i) *= wcum(i);
  }
  arma::mat out(X * arma::trans(sp_out));
  return out;
}


//' Calculating binned matrix product for response term for sparse matrices
//'
//' This function calculates the matrix product (sparse matrices) using Algorithm 3 of Zheyuan Li, Simon N. Wood: "Faster
//' model matrix crossproducts for large generalized linear models with discretized covariates". The idea
//' is to compute just on the unique rows of X by also using an index vector to map to the original matrix.
//' The algorithm implemented here is a small adaption of the original algorithm. Instead of calculating $XW$
//' which again, needs to be transposed, we directly calculate $X^TW$ to avoid another transposing step. In addition
//' to the original algorithm the algorithm here directly calculates the crossproduct with the response.
//'
//' @param X [\code{arma::sp_mat}]\cr
//'   Matrix X.
//' @param y [\code{arma::vec}]\cr
//'   Response vector y.
//' @param k [\code{arma::uvec}]\cr
//'   Index vector for mapping to original matrix $X_o(i,) = X(k(i),.)$.
//' @param w [\code{arma::vec}]\cr
//'   Vector of weights that are accumulated.
//' @return \code{arma::mat} Matrix Product $X^TWX$.
//' @examples
//' nsim = 1e6L
//' nunique = trunc(sqrt(nsim))
//'
//' xunique = runif(n = nunique, min = 0, max = 10)
//' k = sample(x = seq_len(nunique), size = nsim, replace = TRUE)
//'
//' X = poly(x = xunique, degree = 20L)
//' y = runif(nsim)
//'
//' binnedMatMultResponse(X = X, y = y, k = k-1, w = 1)
//' @export
// [[Rcpp::export]]
arma::mat binnedSparseMatMultResponse (const arma::sp_mat& X, const arma::vec& y,  const arma::uvec& k, const arma::vec& w)
{
  unsigned int n = k.size();
  unsigned int ind;

  arma::colvec wcum(X.n_cols, arma::fill::zeros);

  if ( (w.size() == 1) && (w(0) == 1) ) {
    // std::cout << "Weight of size 1:" << std::endl;
    for (unsigned int i = 0; i < n; i++) {
       ind = k(i);
       wcum(ind) += y(i);
    }
  } else {
    for (unsigned int i = 0; i < n; i++) {
      ind = k(i);
      wcum(ind) += w(i) * y(i);
    }
  }
  return X * wcum;
}

#endif // BINNED_MAT_MULT_
