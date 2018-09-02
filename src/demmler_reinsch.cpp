#ifndef DEMMLER_REINSCH_CPP_
#define DEMMLER_REINSCH_CPP_

// [[Rcpp::depends(RcppArmadillo)]]    
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]    
#include <boost/math/tools/toms748_solve.hpp>

#include <functional> // functional

// Define function to calculate the degrees of freedom depending on the penalty term. The root of this function
// corresponds to the desired penalty term:
double calculateDegreesOfFreedom (const double& x, const arma::vec& singular_values, const double& degrees_of_freedom)
{
  return 2 * arma::accu(1 / (1 + x * singular_values)) - arma::accu(1 / arma::pow(1 + x * singular_values, 2)) - degrees_of_freedom;
}

// Use TOMS Algorithm 748: it uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of f(x)
// Included from the boost library: https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots/roots_noderiv/TOMS748.html
double findLambdaWithToms748 (const arma::vec& singular_values, const double& degrees_of_freedom, const double& lower_bound = 0.,
  const double& upper_bound = 1e15) 
{  
  boost::uintmax_t max_iter = 500;
  boost::math::tools::eps_tolerance<double> tol(30);
 
  // Conduct the root finding:
  std::pair<double, double> r = boost::math::tools::toms748_solve(std::bind(calculateDegreesOfFreedom, std::placeholders::_1, singular_values, degrees_of_freedom), lower_bound, upper_bound, tol, max_iter);

  return (r.first + r.second) / 2;
}


//' Transform degrees of freedom to lambda
//'
//' This function calculates the Demmler-Reinsch-Orthogonalization to translate
//' the degrees of freedom to a penalty term.
//'
//' @param XtX [\code{matrix}]\cr
//'   Square matrix calculated by $X^TW^TWX$, where $X$ is the design matrix, and
//'   $W$ the diagonal matrix where the diagonal includes the weights.
//' @param penalty_mat [\code{matrix}]\cr
//'   Penalization matrix used for model fitting.
//' @param degrees_of_freedom [\code{numeric(1)}]\cr
//'   Degrees of freedom to convert to the penalty term.
//' @return \code{numeric(1)} Penalty term which corresponds to the given degrees of freedom.
//' @examples
//' X = cbind(1, iris$Petal.Length, iris$Sepal.Length)
//' weights = rep(1, nrow(iris))
//' pen = penaltyMat(ncol(X), 2)
//' XtX = t(X) %*% t(diag(weights)) %*% diag(weights) %*% X
//' 
//' demmlerReinsch(XtX, pen, 2)
//' @export
// [[Rcpp::export]]
double demmlerReinsch (const arma::mat& XtX, const arma::mat& penalty_mat, const double& degrees_of_freedom)
{
  const double eps = 1e-9;
  // const int XtX_rank = arma::rank(XtX);

  // if (df > XtX_rank) {
  //   Rcpp::stop(\"Degrees of freedom has to be smaller than the rank of the design matrix.\");
  // }
  // if (df == XtX_rank) {
  //   Rcpp::warning(\"Degrees of freedom matches rank of matrix, hence lambda is set to 0.\");
  // }
  arma::mat cholesky = arma::chol(XtX + penalty_mat * eps);
  arma::mat cholesky_inv = arma::inv(cholesky);

  arma::mat Ld  = cholesky_inv.t() * penalty_mat * cholesky_inv;

  arma::vec singular_values = svd(Ld);

  return findLambdaWithToms748(singular_values, degrees_of_freedom);
}

# endif // DEMMLER_REINSCH_CPP_