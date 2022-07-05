#ifndef CENTER_MAT_CPP_
#define CENTER_MAT_CPP_

#include <RcppArmadillo.h>

//' Subtract the basis of one matrix of another
//'
//' This function subtracts the basis of X2 from X1. The resulting
//' matrix is not able to reconstruct the basis of X2.
//' @param X1 [`matrix()`]\cr
//'   Matrix of rank p1 from which the X1 basis is subtracted.
//' @param X2 [`matrix()`]\cr
//'   Matrix of rank p2.
//' @return `matrix()` Matrix of rank p1-p2 with new basis
//' @examples
//' x = runif(100)
//' X1 = cbind(1, x, x^2, x^3, x^4)
//' X2 = cbind(1, x)
//' getSubtractionRotation(X1, X2)
//' @export
// [[Rcpp::export]]
arma::mat getSubtractionRotation (const arma::mat& X1, const arma::mat& X2)
{
  // Cross Product X1 and X2
  arma::mat cross = X1.t() * X2;
  // QR decomp
  // We require and orthogonal matrix Q
  arma::mat R;
  arma::mat Q;
  arma::qr(Q,R,cross);
  // get rank of R and add 1
  int rankR = arma::rank(R);
  // construct Z from rows 0 to last row and column R+1 to last column
  arma::mat Z = Q.cols(rankR,Q.n_cols-1);
  return Z;
}

# endif // CENTER_MAT_CPP_
