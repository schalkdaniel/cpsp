context("Demmler-Reinsch-Orthogonalization works")

test_that("Demmler-Reinsch-Orthogonalization can be computed", {

  X = cbind(1, iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
  D = expect_silent(penaltyMat(ncol(X), 2))

  # Value is obtained from `unname(mboost:::df2lambda(X = X, df = 2.5, dmat = D, weights = 1)["lambda"])`
  df_from_mboost = round(19.703878585735374429, 15)
  expect_equal(round(demmlerReinsch(t(X) %*% X, D, 2.5), 15), df_from_mboost)

  # Speed comparison:
  # microbenchmark::microbenchmark(
  #   "cpsp" = demmlerReinsch(t(X) %*% X, D, 2.5),
  #   "mboost" = mboost:::df2lambda(X = X, df = 2.5, dmat = D, weights = 1)
  # )
  # > Unit: microseconds
  # >    expr    min     lq   mean median     uq    max neval cld
  # >    cpsp  101.4  123.5  164.6  172.1  200.7  287.2   100  a
  # >  mboost 1998.5 2095.4 2269.4 2196.3 2376.3 2867.8   100   b
})
