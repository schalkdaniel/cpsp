context("Tensor products are working")

test_that("Tensors can be calculated", {

  nsim = 100L
  x1 = rnorm(nsim)
  x2 = rnorm(nsim)

  # ------------------------------------------- #
  # Create basis to run algorithms:
  knots1 = createKnots(values = x1, n_knots = 20, degree = 3)
  knots2 = createKnots(values = x2, n_knots = 20, degree = 3)
  X1 = createSplineBasis(values = x1, degree = 3, knots = knots1)
  X2 = createSplineBasis(values = x2, degree = 3, knots = knots2)
  X1_sp = createSparseSplineBasis(values = x1, degree = 3, knots = knots1)
  X2_sp = createSparseSplineBasis(values = x2, degree = 3, knots = knots2)
  # ------------------------------------------- #

  tensor = expect_silent(rowWiseTensor(X1, X2))
  tensor_sp = expect_silent(rowWiseTensorSparse(X1_sp, X2_sp))
})
