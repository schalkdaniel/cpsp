
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R-CMD-check](https://github.com/schalkdaniel/compboostSplines/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/schalkdaniel/compboostSplines/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/schalkdaniel/compboostSplines/branch/main/graph/badge.svg?token=KGR22VAOHI)](https://codecov.io/gh/schalkdaniel/compboostSplines)
[![License: LGPL
v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

## C++ Spline Implementation of Compboost

This repository contains the spline implementation of
[compboost](https://compboost.org) by providing the spline functions
directly without any bloated code around it.

**Feel free to extend the algorithms, improve performance, or use for
your own projets.**

## Installation

#### Developer version:

``` r
devtools::install_github("schalkdaniel/compboostSplines")
```

## Examples

  - [Spline Regression](#spline-regression)
  - [Demmler-Reinsch-Orthogonalization](#demmler-reinsch-orthogonalization)

### Spline Regression

This package can build spline bases for you, ether as dense or sparse
matrix. With the matrix it is possible to do e.g. spline regression or
other cool stuff:

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
basis[1:10, 1:10]
#>              [,1]       [,2]      [,3]         [,4]         [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,] 0.166666667 0.66666667 0.1666667 0.000000e+00 0.0000000000    0    0    0    0     0
#>  [2,] 0.137397362 0.66290097 0.1996613 4.038753e-05 0.0000000000    0    0    0    0     0
#>  [3,] 0.132957646 0.66159342 0.2053853 6.365753e-05 0.0000000000    0    0    0    0     0
#>  [4,] 0.028376319 0.51225485 0.4446072 1.476165e-02 0.0000000000    0    0    0    0     0
#>  [5,] 0.008726429 0.39751898 0.5528898 4.086475e-02 0.0000000000    0    0    0    0     0
#>  [6,] 0.000000000 0.11215656 0.6523140 2.352140e-01 0.0003153798    0    0    0    0     0
#>  [7,] 0.000000000 0.10369088 0.6468246 2.489625e-01 0.0005220619    0    0    0    0     0
#>  [8,] 0.000000000 0.09826255 0.6426961 2.583396e-01 0.0007017948    0    0    0    0     0
#>  [9,] 0.000000000 0.04684948 0.5682074 3.781031e-01 0.0068399859    0    0    0    0     0
#> [10,] 0.000000000 0.03956349 0.5492555 4.019762e-01 0.0092048340    0    0    0    0     0

# You can also create sparse matrices:
basis.sparse = createSparseSplineBasis(values = x, degree = 3, knots = knots)
str(basis.sparse)
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:398] 0 1 2 3 4 0 1 2 3 4 ...
#>   ..@ p       : int [1:25] 0 5 17 33 51 70 82 93 105 115 ...
#>   ..@ Dim     : int [1:2] 100 24
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:398] 0.16667 0.1374 0.13296 0.02838 0.00873 ...
#>   ..@ factors : list()

# Check if row sums add up to 1:
rowSums(basis)
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [48] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [95] 1 1 1 1 1 1

# Polynomial regression using b-splines:
beta = solve(t(basis) %*% basis) %*% t(basis) %*% y

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
beta_pen = solve(t(basis) %*% basis + penalty * K) %*% t(basis) %*% y

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

![](Readme_files/unnamed-chunk-3-1.png)<!-- -->

### Demmler-Reinsch-Orthogonalization

In order to compare different models such as a linear model and additive
model (using splines) we need to set the degrees of freedom equally. The
Demmler-Reinsch-Orthogonalization can be used to translate given degrees
of freedom to a penalty term:

``` r
# We use the basis and penalty matrix from above and specify 2 and 4 degrees of freedom:
(penalty_df2 = demmlerReinsch(t(basis) %*% basis, K, 2))
#> [1] 23719580783
(penalty_df4 = demmlerReinsch(t(basis) %*% basis, K, 4))
#> [1] 390.0661

# This is now used for a new estimator:
beta_df2 = solve(t(basis) %*% basis + penalty_df2 * K) %*% t(basis) %*% y
beta_df4 = solve(t(basis) %*% basis + penalty_df4 * K) %*% t(basis) %*% y

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

![](Readme_files/unnamed-chunk-4-1.png)<!-- -->
