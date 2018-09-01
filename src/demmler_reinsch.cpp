#ifndef DEMMLER_REINSCH_CPP_
#define DEMMLER_REINSCH_CPP_

#include <RcppArmadillo.h>

/* Notes:
 * - Use armadillo svd with method divide-and-conquer to speed up the svd
 * - Which method for root finding? Brent or Newton?
 * - Cholesky decomposition of XtX or XtX + D???
 */

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
//' pen = penaltyMat(10, 2)
//' XtX = t(X) %*% t(diag(weights)) %*% diag(weights) %*% X
//' 
//' demmlerReinsch(XtX, pen, 2)
//' @export
// [[Rcpp::export]]
double demmlerReinsch (const arma::mat& XtX, const arma::mat& penalty_mat, const double& degrees_of_freedom)
{
  // Calculate cholesky decomposition:
  arma::mat cholesky = arma::chol(XtX + penalty_mat)

  // Calculate the difference matrix for higher orders:
  if (differences > 1) {
    arma::mat diffs_reduced = diffs;
    for (unsigned int k = 0; k < differences - 1; k++) {
      diffs_reduced = diffs_reduced(arma::span(1, diffs_reduced.n_rows - 1), arma::span(1, diffs_reduced.n_cols - 1));
      diffs = diffs_reduced * diffs;
    }
  }
  return diffs.t() * diffs;
}



# endif // DEMMLER_REINSCH_CPP_