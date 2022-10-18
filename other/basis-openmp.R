src = '
#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(findSpan)]]
unsigned int findSpan (const double& x, const arma::vec& knots)
{
  // Special case which the algorithm cant handle:
  if (x < knots[1]) { return 0; }
  if (x == knots[knots.size() - 1]) { return knots.size() - 1; }

  unsigned int low = 0;
  unsigned int high = knots.size() - 1;
  unsigned int mid = std::round( (low + high) / 2 );

  while (x < knots[mid] || x >= knots[mid + 1]) {
    if (x < knots[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = std::round( (low + high) / 2 );
  }
  return mid;
}

// [[Rcpp::export(createSplineBasisOMP)]]
arma::mat createSplineBasisOMP (const arma::vec& values, const unsigned int& degree,
  const arma::vec& knots, const int ncores)
{
  unsigned int n_cols =  knots.size() - (degree + 1);
  double x;
  unsigned int idx;

  // Frame for output:
  arma::mat spline_basis(values.size(), n_cols, arma::fill::zeros);

  // Inserting rowwise. This loop creates the basis functions for each row:
  omp_set_num_threads(ncores);
  # pragma omp parallel shared(spline_basis)
  {
    # pragma omp for
    for (unsigned int actual_row = 0; actual_row < values.size(); actual_row++) {

      x = values(actual_row);

      // Index of x within the konts:
      idx = findSpan(x, knots);

      // A problem occurs if x = max(knots), then idx is bigger than
      // number of columns which couses problems. Catch that:
      if (idx > (n_cols - 1)) { idx = n_cols - 1; }

      // Output for basis functions. Here we have the non-zero entries:
      arma::rowvec N(degree + 1, arma::fill::zeros);
      N[0] = 1.0;

      arma::vec left(degree + 1, arma::fill::zeros);
      arma::vec right(degree + 1, arma::fill::zeros);

      double saved;
      double temp;

      // De Boors algorithm to recursive find base in a triangle scheme:
      for (unsigned int j = 1; j <= degree; j++) {

        left[j]  = x - knots[idx + 1 - j];
        right[j] = knots[idx + j] - x;

        saved = 0;

        for (unsigned int r = 0; r < j; r++) {
          temp  = N[r] / (right[r + 1] + left[j - r]);
          N[r]  = saved + right[r + 1] * temp;
          saved = left[j - r] * temp;
        }
        N[j] = saved;
      }
      spline_basis(actual_row, arma::span(idx - degree, idx)) = N;
    }
  }
  return spline_basis;
}


// [[Rcpp::export(createSplineBasis)]]
arma::mat createSplineBasis (const arma::vec& values, const unsigned int& degree,
  const arma::vec& knots)
{
  unsigned int n_cols =  knots.size() - (degree + 1);

  // Index for binary search:
  unsigned int idx;
  // Variable for value on which the basis should be computed:
  double x;

  // Frame for output:
  arma::mat spline_basis(values.size(), n_cols, arma::fill::zeros);

  // Inserting rowwise. This loop creates the basis functions for each row:
  for (unsigned int actual_row = 0; actual_row < values.size(); actual_row++) {

    x = values(actual_row);

    // Index of x within the konts:
    idx = findSpan(x, knots);

    // A problem occurs if x = max(knots), then idx is bigger than
    // number of columns which couses problems. Catch that:
    if (idx > (n_cols - 1)) { idx = n_cols - 1; }

    // Output for basis functions. Here we have the non-zero entries:
    arma::rowvec N(degree + 1, arma::fill::zeros);
    N[0] = 1.0;

    arma::vec left(degree + 1, arma::fill::zeros);
    arma::vec right(degree + 1, arma::fill::zeros);

    double saved;
    double temp;

    // De Boors algorithm to recursive find base in a triangle scheme:
    for (unsigned int j = 1; j <= degree; j++) {

      left[j]  = x - knots[idx + 1 - j];
      right[j] = knots[idx + j] - x;

      saved = 0;

      for (unsigned int r = 0; r < j; r++) {
        temp  = N[r] / (right[r + 1] + left[j - r]);
        N[r]  = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      N[j] = saved;
    }
    spline_basis(actual_row, arma::span(idx - degree, idx)) = N;
  }
  return spline_basis;
}

#include <cmath>

// [[Rcpp::export(getIndex)]]
std::vector<std::vector<unsigned int>> getIndex(unsigned int nc, unsigned int ncores) {
  unsigned int n = nc - 1;
  unsigned int n0 = std::ceil(n / (double)ncores);

  std::vector<std::vector<unsigned int>> out;
  unsigned int upper;
  for (unsigned int i = 0; i < ncores; i++) {
    std::vector<unsigned int> in;
    in.push_back(i * n0);
    upper = (i + 1) * n0 - 1;
    if (upper > n) {
      in.push_back(n);
    } else {
      if (i == (ncores - 1)) {
        in.push_back(n);
      } else {
        in.push_back(upper);
      }
    }
    out.push_back(in);
  }
  return out;
}
// [[Rcpp::export(fillTest)]]
std::vector<std::vector<int>> fillTest(unsigned int n, unsigned int ncores) {

  auto idx = getIndex(n, ncores);

  std::vector<std::vector<int>> out;
  #pragma omp parallel num_threads(ncores) shared(idx)
  {
    int kk = omp_get_thread_num();
    std::vector<int> temp;
    # pragma omp parallel for schedule(static)
    for (unsigned int i = idx[kk][0]; i <= idx[kk][1]; i++) {
      temp.push_back(i);
    }

    # pragma omp critical
    {
      out[kk] = temp;
    }
  }
  return(out);
}
/*
using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export(allFiniteOMP)]]
bool allFiniteOMP(NumericVector x, int ncores)
{
 size_t n = x.size();
 double out = 0;
 NumericVector z(ncores);

 #pragma omp parallel num_threads(ncores)
 {
  int kk = omp_get_thread_num();

  #pragma omp for schedule(static)
  for(size_t ii = 0; ii < n; ii++)
  {
   z[kk] += x[ii];
  }
 }

 out = sum(z);

 return R_FINITE(out);

 }
*/
'


Rcpp::sourceCpp(code = src, verbose = TRUE, rebuild = TRUE)

getIndex(10, 3)

nsim = 1000000L
microbenchmark::microbenchmark(
  test1 = fillTest(nsim, 1),
  test2 = fillTest(nsim, 2),
  test4 = fillTest(nsim, 4),
  times = 3L
)

fillTest(10, 4)

x = runif(1000000, 0, 10)
knots = cpsp::createKnots(x, 20, 3)

#X1 = createSplineBasis(x, 3, knots)
#X2 = createSplineBasisOMP(x, 3, knots, 4)

#all(X1 == X2)

microbenchmark::microbenchmark(
  "sequential" = {createSplineBasis(x, 3, knots)},
  "omp" = {createSplineBasisOMP(x, 3, knots, 4)},
  times = 10L
)
