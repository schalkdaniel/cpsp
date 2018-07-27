[![Build Status](https://travis-ci.org/schalkdaniel/compboostSplines.svg?branch=master)](https://travis-ci.org/schalkdaniel/compboostSplines)
[![Coverage Status](https://coveralls.io/repos/github/schalkdaniel/compboostSplines/badge.svg)](https://coveralls.io/github/schalkdaniel/compboostSplines)


## C++ Spline Implementation of Compboost

This repository contains the spline implementation of compboost by providing the spline functions directly without any bloated code around it. 

**Feel free to extend the algorithms, improve performance, or use for your own projets.**

## Installation

#### Developer version:

```r
devtools::install_github("schalkdaniel/compboostSplines")
```

## Examples

This package can build spline bases for you, ether as dense or sparse matrix. With the matrix it is possible to do e.g. spline regression or other cool stuff:
```r
library(compboostSpliens)

nsim = 100

# Sample data:
x = sort(runif(nsim, 0, 10))
y = 2 * sin(x) + rnorm(nsim, 0, 0.5)

# Calculate knots of given x values:
knots = createKnots(values = x, n_knots = 20, degree = 3)

# Create basis using that knots:
basis = createSplineBasis(values = x, degree = 3, knots = knots)
basis[1:10, 1:10]
##            [,1]     [,2]    [,3]     [,4]    [,5]     [,6] [,7] [,8] [,9] [,10]
##  [1,] 0.1666667 0.666667 0.16667 0.000000 0.00000 0.000000    0    0    0     0
##  [2,] 0.0508134 0.577242 0.36612 0.005825 0.00000 0.000000    0    0    0     0
##  [3,] 0.0435305 0.559983 0.38866 0.007827 0.00000 0.000000    0    0    0     0
##  [4,] 0.0016032 0.290810 0.62625 0.081341 0.00000 0.000000    0    0    0     0
##  [5,] 0.0012435 0.279729 0.63221 0.086813 0.00000 0.000000    0    0    0     0
##  [6,] 0.0002767 0.232053 0.65348 0.114195 0.00000 0.000000    0    0    0     0
##  [7,] 0.0000000 0.016360 0.45463 0.502954 0.02606 0.000000    0    0    0     0
##  [8,] 0.0000000 0.005581 0.36303 0.579521 0.05187 0.000000    0    0    0     0
##  [9,] 0.0000000 0.000000 0.04624 0.566734 0.38002 0.007012    0    0    0     0
## [10,] 0.0000000 0.000000 0.01405 0.439872 0.51657 0.029512    0    0    0     0

# You can also create sparse matrices:
basis.sparse = createSparseSplineBasis(values = x, degree = 3, knots = knots)
str(basis.sparse)
## Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
##   ..@ i       : int [1:398] 0 1 2 3 4 5 6 0 1 2 ...
##   ..@ p       : int [1:25] 0 7 18 37 58 77 101 119 135 153 ...
##   ..@ Dim     : int [1:2] 100 24
##   ..@ Dimnames:List of 2
##   .. ..$ : NULL
##   .. ..$ : NULL
##   ..@ x       : num [1:398] 0.1667 0.1315 0.0691 0.0495 0.0267 ...
##   ..@ factors : list()


# Check if row sums add up to 1:
rowSums(basis)
##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
## [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
## [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1


# Polynomial regression using b-splines:
beta = solve(t(basis) %*% basis) %*% t(basis) %*% y

# 20 knots may tend to overfit on the data, lets try p-splines with a penalty term of 4!
penalty = 4

# Get penalty matrix:
K = penaltyMat(ncol(basis), differences = 2)
K[1:6, 1:6]
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    1   -2    1    0    0    0
## [2,]   -2    5   -4    1    0    0
## [3,]    1   -4    6   -4    1    0
## [4,]    0    1   -4    6   -4    1
## [5,]    0    0    1   -4    6   -4
## [6,]    0    0    0    1   -4    6


# Get new estimator:
beta.pen = solve(t(basis) %*% basis + penalty * K) %*% t(basis) %*% y

# Lets visualize the curves:

library(ggplot2)
library(ggthemes)

plot.df = data.frame(
	x = x,
	y = y
)
spline.df = data.frame(
	"Spline" = c(rep("B-Spline", nsim), rep("P-Spline", nsim)),
	"x" = rep(x, 2),
	"y" = c(basis %*% beta, basis %*% beta.pen)
)

ggplot() + geom_point(data = plot.df, mapping = aes(x = x, y = y)) +
  geom_line(data = spline.df, mapping = aes(x = x, y = y, color = Spline)) +
  theme_tufte() + 
  scale_color_brewer(palette = "Set1")
```
<p align="center">
  <img src="other/spline.png?raw=true" alt="Spline Visualization" width="70%">
</p>


## License

© 2018 [Daniel Schalk](https://danielschalk.com)

The contents of this repository are distributed under the MIT license. See below for details:

> The MIT License (MIT)
> 
> Copyright (c) 2018 Daniel Schalk
> 
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
> 
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
> 
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
