src1 = "
double myTest1 (arma::vec x, double a)
{
  return arma::accu(1 / arma::pow(1 + a * x, 2));
}
"

src2 = "
double myTest2 (const double& lambda, const arma::vec& singular_values) 
{
  return 2 * arma::accu(1 / (1 + lambda * singular_values)) - arma::accu(1 / arma::pow(1 + lambda * singular_values, 2));
}
"

src_demmler_reinsch = "
arma::vec demmlerReinsch (const arma::mat& XtX, const arma::mat& penalty_mat, const double& degrees_of_freedom)
{
  const double eps = 1e-9;
  const int XtX_rank = arma::rank(XtX);

  // if (df > XtX_rank) {
  //   Rcpp::stop(\"Degrees of freedom has to be smaller than the rank of the design matrix.\");
  // }
  // if (df == XtX_rank) {
  //   Rcpp::warning(\"Degrees of freedom matches rank of matrix, hence lambda is set to 0.\");
  // }
  // Calculate cholesky decomposition:
  arma::mat cholesky = arma::chol(XtX + penalty_mat * eps);
  arma::mat cholesky_inv = arma::inv(cholesky);

  arma::mat Ld  = cholesky_inv.t() * penalty_mat * cholesky_inv;

  arma::vec singular_values = svd(Ld);

  return singular_values;
}
"

Rcpp::cppFunction(code = src1, depends = "RcppArmadillo")
Rcpp::cppFunction(code = src2, depends = "RcppArmadillo")
# Rcpp::cppFunction(code = src_demmler_reinsch, depends = "RcppArmadillo")

# demmlerReinsch(XtX, D, 2)


Rcpp::sourceCpp("other/test_brent.cpp")

myFun(2, 1:3, 4)

x_plot = seq(-2, 5, 0.01)
y = exp(x_plot) - x_plot^2

plot(x = x_plot, y = y, type = "l")

optimBrent(1:3, 4, lower_bound = 0, upper_bound = 1e15)



