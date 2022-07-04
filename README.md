
<!-- README.md is generated from README.Rmd. Please edit that file -->

    ## Rcpp         (NA -> 1.0.8.3   ) [CRAN]
    ## BH           (NA -> 1.78.0-0  ) [CRAN]
    ## RcppArmad... (NA -> 0.11.2.0.0) [CRAN]

    ## Installing 3 packages: Rcpp, BH, RcppArmadillo

    ## Installing packages into '/Users/runner/work/_temp/Library'
    ## (as 'lib' is unspecified)

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T//RtmpzCHz4y/downloaded_packages
    ## * checking for file ‘/Users/runner/work/compboostSplines/compboostSplines/DESCRIPTION’ ... OK
    ## * preparing ‘compboostSplines’:
    ## * checking DESCRIPTION meta-information ... OK
    ## * cleaning src
    ## * checking for LF line-endings in source and make files and shell scripts
    ## * checking for empty or unneeded directories
    ## Omitted ‘LazyData’ from DESCRIPTION
    ## * building ‘compboostSplines_0.1.tar.gz’
    ## 
    ## Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
    ##   /var/folders/24/8k48jl6d249_n_qfxwsl6xvm0000gn/T//RtmpzCHz4y/compboostSplines_0.1.tar.gz \
    ##   --install-tests 
    ## * installing to library ‘/Users/runner/work/_temp/Library’
    ## * installing *source* package ‘compboostSplines’ ...
    ## ** using staged installation
    ## ** libs
    ## clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/runner/work/_temp/Library/Rcpp/include' -I'/Users/runner/work/_temp/Library/RcppArmadillo/include' -I'/Users/runner/work/_temp/Library/BH/include' -I/usr/local/include    -fPIC  -Wall -g -O2  -c RcppExports.cpp -o RcppExports.o
    ## clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/runner/work/_temp/Library/Rcpp/include' -I'/Users/runner/work/_temp/Library/RcppArmadillo/include' -I'/Users/runner/work/_temp/Library/BH/include' -I/usr/local/include    -fPIC  -Wall -g -O2  -c center_matrices.cpp -o center_matrices.o
    ## clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/runner/work/_temp/Library/Rcpp/include' -I'/Users/runner/work/_temp/Library/RcppArmadillo/include' -I'/Users/runner/work/_temp/Library/BH/include' -I/usr/local/include    -fPIC  -Wall -g -O2  -c demmler_reinsch.cpp -o demmler_reinsch.o
    ## clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/runner/work/_temp/Library/Rcpp/include' -I'/Users/runner/work/_temp/Library/RcppArmadillo/include' -I'/Users/runner/work/_temp/Library/BH/include' -I/usr/local/include    -fPIC  -Wall -g -O2  -c splines.cpp -o splines.o
    ## clang++ -mmacosx-version-min=10.13 -std=gnu++11 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Users/runner/work/_temp/Library/Rcpp/include' -I'/Users/runner/work/_temp/Library/RcppArmadillo/include' -I'/Users/runner/work/_temp/Library/BH/include' -I/usr/local/include    -fPIC  -Wall -g -O2  -c tensor.cpp -o tensor.o
    ## clang++ -mmacosx-version-min=10.13 -std=gnu++11 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o compboostSplines.so RcppExports.o center_matrices.o demmler_reinsch.o splines.o tensor.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin18/8.2.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
    ## ld: warning: dylib (/usr/local/gfortran/lib/libgfortran.dylib) was built for newer macOS version (10.14) than being linked (10.13)
    ## ld: warning: dylib (/usr/local/gfortran/lib/libquadmath.dylib) was built for newer macOS version (10.14) than being linked (10.13)
    ## installing to /Users/runner/work/_temp/Library/00LOCK-compboostSplines/00new/compboostSplines/libs
    ## ** R
    ## ** tests
    ## ** byte-compile and prepare package for lazy loading
    ## ** help
    ## *** installing help indices
    ## ** building package indices
    ## ** testing if installed package can be loaded from temporary location
    ## ** checking absolute paths in shared objects and dynamic libraries
    ## ** testing if installed package can be loaded from final location
    ## ** testing if installed package keeps a record of temporary installation path
    ## * DONE (compboostSplines)

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
#>               [,1]       [,2]      [,3]       [,4]         [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,] 1.666667e-01 0.66666667 0.1666667 0.00000000 0.0000000000    0    0    0    0     0
#>  [2,] 6.692429e-02 0.60691150 0.3231583 0.00300590 0.0000000000    0    0    0    0     0
#>  [3,] 1.203597e-02 0.42548261 0.5293589 0.03312249 0.0000000000    0    0    0    0     0
#>  [4,] 1.180423e-03 0.27758648 0.6333278 0.08790530 0.0000000000    0    0    0    0     0
#>  [5,] 2.011503e-04 0.22496566 0.6559344 0.11889876 0.0000000000    0    0    0    0     0
#>  [6,] 7.395761e-05 0.20749081 0.6610708 0.13136446 0.0000000000    0    0    0    0     0
#>  [7,] 6.314715e-05 0.20527551 0.6616201 0.13304128 0.0000000000    0    0    0    0     0
#>  [8,] 0.000000e+00 0.12770600 0.6597598 0.21243208 0.0001021025    0    0    0    0     0
#>  [9,] 0.000000e+00 0.06267074 0.6000287 0.33371129 0.0035893005    0    0    0    0     0
#> [10,] 0.000000e+00 0.02679960 0.5060105 0.45136428 0.0158255883    0    0    0    0     0

# You can also create sparse matrices:
basis.sparse = createSparseSplineBasis(values = x, degree = 3, knots = knots)
str(basis.sparse)
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   ..@ i       : int [1:398] 0 1 2 3 4 5 6 0 1 2 ...
#>   ..@ p       : int [1:25] 0 7 17 27 40 51 65 84 106 128 ...
#>   ..@ Dim     : int [1:2] 100 24
#>   ..@ Dimnames:List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   ..@ x       : num [1:398] 0.166667 0.066924 0.012036 0.00118 0.000201 ...
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
#> [1] 46980601292
(penalty_df4 = demmlerReinsch(t(basis) %*% basis, K, 4))
#> [1] 448.0386

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
