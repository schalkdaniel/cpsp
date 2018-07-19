context("splines can be computed")

test_that("splines can be computed", {

  nsim = 100000

  x = sort(runif(nsim, 0, 10))
  y = 2 * sin(x) + rnorm(nsim, 0, 0.5)

  expect_silent(penaltyMat(10, 2))

  expect_silent({
    knots = createKnots(values = x, n_knots = 100, degree = 3)
  })
  expect_silent({
    X = createSplineBasis(values = x, degree = 3, knots = knots)
  })
  expect_silent({
    X.sparse = createSparseSplineBasis(values = x, degree = 3, knots = knots)
  })

  X.sparse.dense = as.matrix(X.sparse)
  attributes(X.sparse.dense) = NULL
  dim(X.sparse.dense) = dim(X)

  expect_equal(X, X.sparse.dense)
})