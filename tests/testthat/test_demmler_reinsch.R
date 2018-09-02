context("Demmler-Reinsch-Orthogonalization works")

test_that("Demmler-Reinsch-Orthogonalization can be computed", {

  X = cbind(1, iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
  D = compboostSplines::penaltyMat(ncol(X), 2)
  
  expect_equal(demmlerReinsch(t(X) %*% X, D, 2.5), unname(mboost:::df2lambda(X = X, df = 2.5, dmat = D, weights = 1)["lambda"]))
})