context("Binning works")

test_that("Binning functions can be run", {

  nsim = 10000L
  x = rnorm(nsim)

  bins1 = expect_silent(binVectorCustom(x, 100))
  bins2 = expect_silent(binVector(x))
  expect_equal(bins1, bins2)
  expect_equal(length(bins1), 100)

  idx = expect_silent(calculateIndexVector(x, bins1))
  expect_equal(min(idx), 0)
  expect_equal(max(idx), 99)

  # ------------------------------------------- #
  # Create basis to run algorithms:
  knots = createKnots(values = x, n_knots = 20, degree = 3)
  basis_bin = createSplineBasis(values = bins1, degree = 3, knots = knots)
  basis_bin_sp = createSparseSplineBasis(values = bins1, degree = 3, knots = knots)
  y = rnorm(nsim)
  w = rep(1, nsim)
  # ------------------------------------------- #

  xtx1 = expect_silent(binnedSparseMatMult(Matrix::t(basis_bin_sp), idx, w))
  xty1 = expect_silent(binnedSparseMatMultResponse(Matrix::t(basis_bin_sp), y, idx, w))
  xtx2 = expect_silent(binnedMatMult(basis_bin, idx, w))
  xty2 = expect_silent(binnedMatMultResponse(basis_bin, y, idx, w))

  expect_equal(xtx1, xtx2)
  expect_equal(xty1, t(xty2))

  basis = createSplineBasis(values = x, degree = 3, knots = knots)
  xtx0 = binnedMatMult(basis, seq_len(nsim) - 1, w)
  xty0 = binnedMatMultResponse(basis, y, seq_len(nsim) - 1, w)

  expect_equal(xtx0, t(basis) %*% basis)
  expect_equal(t(xty0), t(basis) %*% y)
})
