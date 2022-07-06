
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R-CMD-check](https://github.com/schalkdaniel/compboostSplines/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/schalkdaniel/compboostSplines/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/schalkdaniel/compboostSplines/branch/main/graph/badge.svg?token=KGR22VAOHI)](https://codecov.io/gh/schalkdaniel/compboostSplines)
[![License: LGPL
v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

## C++ Spline Implementation of Compboost

Pure spline implementation of [compboost](https://compboost.org). This
repository can be used for spline regression and non-parametric
modelling of numerical features. The functionality contains the basis
creation, the Demmler-Reinsch-Orthogonalization, as well as a method to
subtract matrix bases.

**Feel free to extend the algorithms, improve performance, or use for
your own projets.**

## Installation

#### Developer version:

``` r
devtools::install_github("schalkdaniel/compboostSplines")
```

## Examples

<!--
- [Spline basis](#spline-basis)
- [Spline regression](#spline-regression)
- [Demmler-Reinsch-Orthogonalization](#demmler-reinsch-orthogonalization)
-->

### Spline basis

To create a spline basis call `createSplineBasis()` (for dense matrices)
and `createSparseSplineBasis()` (for sparse matrices):

``` r
library(compboostSplines)

nsim = 100

# Sample data:
x = sort(runif(nsim, 0, 10))
y = 2 * sin(x) + rnorm(nsim, 0, 0.5)

# Calculate knots of given x values:
knots = createKnots(values = x, n_knots = 20, degree = 3)

# Create basis using that knots:
basis = createSplineBasis(values = x, degree = 3, knots = knots)
str(basis)
#>  num [1:100, 1:24] 0.16667 0.03642 0.01341 0.00912 0.00118 ...

# You can also create sparse matrices:
basis_sparse = createSparseSplineBasis(values = x, degree = 3, knots = knots)
str(basis_sparse)
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:398] 0 1 2 3 4 0 1 2 3 4 ...
#>   ..@ p       : int [1:25] 0 5 13 24 41 62 85 107 128 145 ...
#>   ..@ Dim     : int [1:2] 100 24
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:398] 0.16667 0.03642 0.01341 0.00912 0.00118 ...
#>   ..@ factors : list()

# Check if row sums add up to 1:
rowSums(basis)
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [48] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [95] 1 1 1 1 1 1
```

### Spline regression

With that, it is straight forward to use the least squares estimator for
spline regression:

``` r
# Custom function to estimate the least squares estimator
myEstimator = function(X, y, penmat = 0) {
  xtx = crossprod(X)
  xty = crossprod(X, y)
  L = chol(xtx + penmat)
  z = backsolve(L, xty, transpose = TRUE)
  return(as.vector(backsolve(L, z)))
}

# Polynomial regression using b-splines:
beta = myEstimator(basis, y)

# 20 knots may tend to overfit on the data, lets try p-splines with a penalty term of 4!
penalty = 4

# Get penalty matrix:
K = penaltyMat(ncol(basis), differences = 2)
K[1:6, 1:6]
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1   -2    1    0    0    0
#> [2,]   -2    5   -4    1    0    0
#> [3,]    1   -4    6   -4    1    0
#> [4,]    0    1   -4    6   -4    1
#> [5,]    0    0    1   -4    6   -4
#> [6,]    0    0    0    1   -4    6

# Get new estimator:
beta_pen = myEstimator(basis, y, penalty * K)

# Lets visualize the curves:

library(ggplot2)
library(ggthemes)

plot_df = data.frame(
  x = x,
  y = y
)
spline_df = data.frame(
  "Spline" = c(rep("B-Spline", nsim), rep("P-Spline", nsim)),
  "x" = rep(x, 2),
  "y" = c(basis %*% beta, basis %*% beta_pen)
)

ggplot() + geom_point(data = plot_df, mapping = aes(x = x, y = y)) +
  geom_line(data = spline_df, mapping = aes(x = x, y = y, color = Spline)) +
  theme_tufte() +
  scale_color_brewer(palette = "Set1")
```

![](Readme_files/unnamed-chunk-4-1.png)<!-- -->

### Demmler-Reinsch-Orthogonalization

In order to compare different models such as a linear model and additive
model (using splines) we need to set the degrees of freedom equally. The
Demmler-Reinsch-Orthogonalization can be used to translate given degrees
of freedom to a penalty term:

``` r
# We use the basis and penalty matrix from above and specify 2 and 4 degrees of freedom:
(penalty_df2 = demmlerReinsch(t(basis) %*% basis, K, 2))
#> [1] 1.051e+11
(penalty_df4 = demmlerReinsch(t(basis) %*% basis, K, 4))
#> [1] 436.7

# This is now used for a new estimator:
beta_df2 = myEstimator(basis, y, penalty_df2 * K)
beta_df4 = myEstimator(basis, y, penalty_df4 * K)

plot_df = data.frame(
  x = x,
  y = y
)
types = c("B-Spline", "P-Spline", "P-Spline with 2 df", "P-Spline with 4 df")
spline_df = data.frame(
  "Spline" = rep(types, each = nsim),
  "x" = rep(x, length(types)),
  "y" = c(basis %*% beta, basis %*% beta_pen, basis %*% beta_df2, basis %*% beta_df4)
)

ggplot() + geom_point(data = plot_df, mapping = aes(x = x, y = y)) +
  geom_line(data = spline_df, mapping = aes(x = x, y = y, color = Spline)) +
  theme_tufte() +
  scale_color_brewer(palette = "Set1")
```

![](Readme_files/unnamed-chunk-5-1.png)<!-- -->

### Subtract model effects

The idea of subtracting effects is to not let a design matrix model the
effect represented in another matrix. For example, let’s define a linear
effect:

``` r
y_linear = 2 + 1.5 * x + rnorm(nsim)
df_linear = data.frame(x = x, y = y_linear)

ggplot() +
  geom_point(data = df_linear, aes(x = x, y = y)) +
  theme_tufte()
```

![](Readme_files/unnamed-chunk-6-1.png)<!-- -->

The linear effect is expressed in the design matrix `(1, x)`. With
`centerMatrix()` we can generate a rotation matrix which is used to
rotate the design matrix that the linear effect cannot modelled anymore
while the full design matrix still models the linear effect:

``` r
X_linear = cbind(1, x)
X_nonlinear = basis %*% getSubtractionRotation(basis, X_linear)

pred_full = basis %*% myEstimator(basis, y_linear)
pred_nonlinear = X_nonlinear %*% myEstimator(X_nonlinear, y_linear)

df_plot = data.frame(
  x = rep(x, 2),
  y = c(pred_nonlinear, pred_full),
  type = rep(c("Non linear effect", "Full effect"), each = length(x)))

ggplot() +
  geom_point(data = data.frame(x = x, y = y_linear), aes(x = x, y = y)) +
  geom_line(data = df_plot, aes(x = x, y = y, color = type)) +
  theme_tufte() +
  scale_color_brewer(palette = "Set1")
```

![](Readme_files/unnamed-chunk-7-1.png)<!-- -->

If we use the non-linear effect from the first examples, we observe that
both, the reduced and full design matrix can capture the non-linear
effect. The only differences observable is the shift on the y axis which
results from the subtraction of a constant, meaning that the estimated
effect is always centered around zero if the design matrix of the linear
effect contains a constant column:

``` r
beta_full_nl = myEstimator(basis, y)
beta_nonlinear_nl = myEstimator(X_nonlinear, y)

pred_full_nl = basis %*% beta_full_nl
pred_nonlinear_nl = X_nonlinear %*% beta_nonlinear_nl

df_plot_nl = data.frame(
  x = rep(x, 2),
  y = c(pred_nonlinear_nl, pred_full_nl),
  type = rep(c("Non linear effect", "Full effect"), each = length(x)))

ggplot() +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_line(data = df_plot_nl, aes(x = x, y = y, color = type)) +
  theme_tufte() +
  scale_color_brewer(palette = "Set1")
```

![](Readme_files/unnamed-chunk-8-1.png)<!-- -->

### Binning

Using binning allows to reduce the number of observations in the design
matrix. Therefore, the vector `x` is discretized into a smaller amount
of bins. Furthermore, a mapping keeps the information, which of the bins
represents which original value. Efficient algorithms are design based
on this information:

``` r
bins = binVectorCustom(x, 50)
idx = calculateIndexVector(x, bins) + 1 # Note the shift from C++ to R in indexing
head(data.frame(x = x, bins = bins[idx]))
#>         x    bins
#> 1 0.02933 0.02933
#> 2 0.21807 0.23274
#> 3 0.29906 0.23274
#> 4 0.32377 0.23274
#> 5 0.41275 0.43615
#> 6 0.51136 0.43615
```

For spline regression, we can build the basis just using the bins and
use this (much) smaller design matrix with the optimized algorithms:

``` r
library(Matrix)

w = rep(1, length(idx))
basis_bin = createSparseSplineBasis(values = bins, degree = 3, knots = knots)
xtx = binnedSparseMatMult(t(basis_bin), idx - 1, w)
beta_bin = solve(xtx) %*% binnedSparseMatMultResponse(t(basis_bin), y, idx - 1, w)

df_bin = data.frame(
  x = rep(x, 2),
  y = c(basis %*% beta, basis %*% beta_bin),
  binning = rep(c("No Binning", "Binning"), each = length(x)))

ggplot() + geom_point(data = plot_df, mapping = aes(x = x, y = y)) +
  geom_line(data = df_bin, mapping = aes(x = x, y = y, color = binning)) +
  theme_tufte() +
  scale_color_brewer(palette = "Set1")
```

![](Readme_files/unnamed-chunk-10-1.png)<!-- -->
