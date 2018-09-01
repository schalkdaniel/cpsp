src = "
arma::mat myChol (arma::mat& A)
{
  return arma::chol(A);
}
"

Rcpp::cppFunction(code = src, depends = "RcppArmadillo")
