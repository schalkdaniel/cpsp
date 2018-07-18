// #ifndef SPLINE_SPARSE_CPP_
// #define SPLINE_SPARSE_CPP_

// #include <RcppArmadillo.h>
// #include "splines.h"

// //' De Boors algorithm to find basis functions
// //' 
// //' @param x `double` Point to search for position in knots.
// //' @param degree `unsigned int` Degree of the polynomial between the knots.
// //' @param knots `arma::vec` Vector of knots. It's the users responsibility to
// //'   pass a **SORTED** vector.
// //'   
// //' @return `arma::sp_mat` Sparse matrix containing the basis functions.
// //' @examples
// //' x = sort(runif(100, 0, 10))
// //' y = 2 * sin(x) + rnorm(100, 0, 0.5)
// //' 
// //' # Create knots on the space of x:
// //' knots = createKnots(values = x, n_knots = 7, degree = 3)
// //' 
// //' # Create basis functions for one value:
// //' sparseBasisFuns(x = x[30], degree = 3, knots = knots)
// //' @export
// // [[Rcpp::export]]
// arma::sp_mat sparseBasisFuns (const double& x, const unsigned int& degree, 
//   const arma::vec& knots)
// {
//   // Sparse output vector of bases:
//   arma::rowvec full_base(knots.size() - (degree + 1), arma::fill::zeros);
  
//   // Index of x within the konts:
//   unsigned int idx = findSpan(x, knots);
  
//   // A problem occurs if x = max(knots), then idx is bigger than
//   // length(full_base) which couses problems. Catch that:
//   if (idx > (full_base.n_cols - 1)) { idx = full_base.n_cols - 1; }
  
//   // Output for basis functions. Here we have the non-zero entries:
//   arma::rowvec N(degree + 1, arma::fill::zeros);
//   N[0] = 1.0;
  
//   arma::vec left(degree + 1, arma::fill::zeros);
//   arma::vec right(degree + 1, arma::fill::zeros);
  
//   double saved;
//   double temp;
  
//   // De Boors algorithm to recursive find base in a triangle scheme:
//   for (unsigned int j = 1; j <= degree; j++) {

//     left[j]  = x - knots[idx + 1 - j];
//     right[j] = knots[idx + j] - x;

//     saved = 0;

//     for (unsigned int r = 0; r < j; r++) {
//       temp  = N[r] / (right[r + 1] + left[j - r]);
//       N[r]  = saved + right[r + 1] * temp;
//       saved = left[j - r] * temp;
//     }
//     N[j] = saved;
//   }
  
//   full_base(arma::span(idx - degree, idx)) = N;
  
//   arma::sp_mat out(1, full_base.size());
//   out.row(0) = full_base.t();
  
//   return out;
// }

// //' Transformation from a vector of input points to sparse matrix of basis
// //' 
// //' This functions takes a vector of points and create a sparse matrix of
// //' basis functions. Each row contains the basis of the corresponding value 
// //' in `values`.
// //' 
// //' @param values `arma::vec` Points to create the basis matrix.
// //' @param n_knots `unsigned int` Number of innter knots.
// //' @param degree `unsigned int` polynomial degree of splines.
// //'    
// //' @return `sp_mat` sparse matrix of base functions.
// //' @examples
// //' nsim = 100
// //' 
// //' x = sort(runif(nsim, 0, 10))
// //' y = 2 * sin(x) + rnorm(nsim, 0, 0.5) 
// //' knots = createKnots(values = x, n_knots = 20, degree = 3)
// //'
// //' # Create spline basis:
// //' basis = createSparseBasis(values = x, degree = 3, knots = knots)
// //' @export
// // [[Rcpp::export]]
// arma::sp_mat createSparseBasis (const arma::vec& values, const unsigned int& degree, 
//   const arma::vec& knots)
// {
//   // Frame for output:
//   arma::sp_mat spline_basis(values.size(), knots.size() - (degree + 1));
  
//   // Fill frame with functions:
//   for (unsigned int i = 0; i < values.size(); i++) {
//     spline_basis.row(i) = sparseBasisFuns(values[i], degree, knots);
//   }
  
//   return spline_basis;
// }

// # endif // SPLINE_SPARSE_CPP_