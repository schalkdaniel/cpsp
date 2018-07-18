nsim = 100

x = sort(runif(nsim, 0, 10))
y = 2 * sin(x) + rnorm(nsim, 0, 0.5)

# Create spline basis:
knots = createKnots(values = x, n_knots = 20, degree = 3)
basisFuns(x = x[30], degree = 3, knots = knots)
basis = createBasis(values = x, degree = 3, knots = knots)

# Check if row sums add up to 1:
rowSums(basis)

# Polynomial regression using b-splines:
beta = solve(t(basis) %*% basis) %*% t(basis) %*% y

plot(x = x, y = y)
points(x = x, y = basis %*% beta, type = "l", col = "red")

# 20 knots seems to much, lets try p-splines with a penalty term of 4!
penalty = 4

# Get penalty matrix:
K = penaltyMat(ncol(basis), differences = 2)

# Get new estimator:
beta.pen = solve(t(basis) %*% basis + penalty * K) %*% t(basis) %*% y

points(x = x, y = basis %*% beta.pen, type = "l", col = "dark green")

# Looks much better now!

# Here a nicer one with ggplot:

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








# Test sparse matrices:

nsim = 10000

x = sort(runif(nsim, 0, 10))
y = 2 * sin(x) + rnorm(nsim, 0, 0.5)

# Create spline basis:
knots = createKnots(values = x, n_knots = 20, degree = 3)

X = createBasis(values = x, degree = 3, knots = knots)
X.sparse = createSparseBasis(values = x, degree = 3, knots = knots)

X.sparse.dense = as.matrix(X.sparse)
attributes(X.sparse.dense) = NULL
dim(X.sparse.dense) = dim(X)

all.equal(X, X.sparse.dense)

(bm = microbenchmark::microbenchmark(
	"dense" = createBasis(values = x, degree = 3, knots = knots),
	"sparse" = sparseBasisFuns(values = x, degree = 3, knots = knots), 
	times = 10L
))



plot(bm)




# sp_mat armaEx(S4 mat, bool show) {
# IntegerVector dims = mat.slot("Dim");
# arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
# arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
# arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));
# int nrow = dims[0], ncol = dims[1];
# arma::sp_mat res(i, p, x, nrow, ncol);
# if (show) Rcpp::Rcout << res << std::endl;
# return res;

test.source = "
arma::sp_mat test (arma::vec i_rows, arma::vec j_cols, arma::vec values, int nrows, int ncols) 
{
	arma::umat idx(2, i_rows.size(), arma::fill::zeros);
	for (int i = 0; i < i_rows.size(); i++) {
		idx(0, i) = i_rows(i);
		idx(1, i) = j_cols(i);
	}
	arma::sp_mat out(idx, values, nrows, ncols);
	return out;
	// return idx;
}
"

Rcpp::cppFunction(code = test.source, depends = "RcppArmadillo")

(X = test(c(2, 0, 1, 2, 0, 1), c(0, 1, 2, 4, 4, 6), c(1, 2, 3, 4, 5, 6), 3, 9))
