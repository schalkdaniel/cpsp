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

  xtx1w = expect_silent(binnedSparseMatMult(Matrix::t(basis_bin_sp), idx, 1))
  xty1w = expect_silent(binnedSparseMatMultResponse(Matrix::t(basis_bin_sp), y, idx, 1))
  xtx2w = expect_silent(binnedMatMult(basis_bin, idx, 1))
  xty2w = expect_silent(binnedMatMultResponse(basis_bin, y, idx, 1))

  xtx1fast = expect_silent(binnedMatMult(basis_bin, idx, w, use_fast_acc = TRUE))
  xtx2fast = expect_silent(binnedMatMult(basis_bin, idx, 1, use_fast_acc = TRUE))

  expect_equal(xtx1, xtx2)
  expect_equal(xtx1, xtx1w)
  expect_equal(xtx1, xtx2w)
  expect_equal(xty1, t(xty2))
  expect_equal(xty1, xty1w)
  expect_equal(xty1, t(xty2w))
  expect_equal(xtx1, xtx1fast)
  expect_equal(xtx1, xtx2fast)

  basis = createSplineBasis(values = x, degree = 3, knots = knots)
  xtx0 = binnedMatMult(basis, seq_len(nsim) - 1, w)
  xty0 = binnedMatMultResponse(basis, y, seq_len(nsim) - 1, w)

  expect_equal(xtx0, t(basis) %*% basis)
  expect_equal(t(xty0), t(basis) %*% y)

  # Microbenchmark for fast accumulattion:
  # microbenchmark::microbenchmark(
  #   "fast" = binnedMatMult(basis_bin, idx, 1, use_fast_acc = TRUE),
  #   "standard" = binnedMatMult(basis_bin, idx, 1)
  # )
  # > Unit: microseconds
  # >      expr    min    lq  mean median     uq   max neval cld
  # >      fast  59.27  66.8  74.8  71.49  84.98 108.0   100  a
  # >  standard 177.28 201.4 225.6 219.21 248.50 355.6   100   b
})
