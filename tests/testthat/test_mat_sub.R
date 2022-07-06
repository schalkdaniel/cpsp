context("Matrix rotation works")

test_that("Matrix rotation can be computed", {

  nsim = 1000L
  x = rnorm(nsim)

  # ------------------------------------------- #
  # Create basis to run algorithms:
  knots = createKnots(values = x, n_knots = 20, degree = 3)
  basis = createSplineBasis(values = x, degree = 3, knots = knots)
  basis_lin = cbind(1, x)
  # ------------------------------------------- #

  rotation = expect_silent(getSubtractionRotation(basis, basis_lin))
  expect_equal(ncol(rotation), ncol(basis) - ncol(basis_lin))
  basis_rot = expect_silent(basis %*% rotation)
  expect_true(all(abs(crossprod(basis_rot, basis_lin)) < 1e-10))
})
