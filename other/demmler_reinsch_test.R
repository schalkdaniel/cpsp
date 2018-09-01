X = cbind(1, iris$Sepal.Length, iris$Sepal.Width)
D = compboostSplines::penaltyMat(3, 2)

debug(mboost:::df2lambda)
mboost:::df2lambda(X = X, df = 2.5, dmat = D, weights = 1)
undebug(mboost:::df2lambda)

df = 2.5
dmat = D
weights = 1

if (df > Matrix::rankMatrix(X)) {
  stop("Degrees of freedom has to be smaller than the rank of the design matrix.")
}
if (df == Matrix::rankMatrix(X)) {
  warning("Degrees of freedom matches rank of matrix, hence lambda is set to 0.")
}

XtX = crossprod(X * sqrt(weights))

A   = XtX + dmat * options("mboost_eps")[[1]]
R   = chol(A)
Ld  = solve(t(R)) %*% D %*% solve(R)
# Rm  = backsolve(chol(A), x = diag(ncol(XtX)))
d   = svd(Ld, nu = 0, nv = 0)$d

check = options("mboost_dftraceS")[[1]]
if (check) {
  dfFun = function(lambda) sum(1/(1 + lambda * d))
} else {
  dfFun = function(lambda) 2 * sum(1/(1 + lambda * d)) - sum(1/(1 + lambda * d)^2)
}

uniroot(f = function (x) { dfFun(x) - df }, interval = c(0, 1e15))

f  = function (x) { dfFun(x) - df }
g  = function (x) { -2 * sum((1 + x * d)^{-2}) + 2 * sum((1 + x * d)^{-3} * d) }

rootNewton = function (fun, grad, x_start, eps = 1e-6, iter_max = 1e3) 
{
  for (k in seq_len(iter_max)) {
    x_new = x_start - fun(x_start) / grad(x_start)
    
    if (abs(x_new - x_start) <= eps) {
      message("Stop root finding by epsilon criteria with a root of ", round(x_new, 4))
      break
    }
    cat("Iteration ", k, ": x = ", round(x_new, 4), ", f(x) = ", round(fun(x_new), 4), ", grad(x) = ", round(grad(x_new)), "\n", sep = "")
    x_start = x_new
  }
  return (x_start)
}

(root = rootNewton(f, g, 3))

x = seq(-1, 4, 0.01)
y = x

for (i in seq_along(x)) {
  y[i] = f(x[i])
}

plot(x = x, y = y, type = "l")
abline(h = 0, col = "red")
abline(v = root, col = "blue")
