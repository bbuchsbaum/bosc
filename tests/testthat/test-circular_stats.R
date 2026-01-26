test_that("circ_mean computes correct mean direction", {
  # All angles at 0 -> mean should be 0
  angles <- rep(0, 10)
  expect_equal(circ_mean(angles), 0, tolerance = 1e-10)

  # All angles at pi/2 -> mean should be pi/2
  angles <- rep(pi / 2, 10)
  expect_equal(circ_mean(angles), pi / 2, tolerance = 1e-10)

  # Symmetric around 0 -> mean should be 0
  angles <- c(-pi / 4, pi / 4)
  expect_equal(circ_mean(angles), 0, tolerance = 1e-10)
})

test_that("circ_mean handles weighted data", {
  angles <- c(0, pi / 2)
  # Equal weights
  expect_equal(circ_mean(angles, w = c(1, 1)), pi / 4, tolerance = 1e-10)
  # Heavier weight on first angle
  result <- circ_mean(angles, w = c(3, 1))
  expect_true(result > 0 && result < pi / 4)
})

test_that("circ_r returns correct resultant length", {
  # Perfectly aligned data -> r = 1
  angles <- rep(0, 10)
  expect_equal(circ_r(angles), 1, tolerance = 1e-10)

  # Uniformly distributed -> r close to 0
  set.seed(42)
  angles_uniform <- seq(0, 2 * pi - 0.01, length.out = 100)
  r <- circ_r(angles_uniform)
  expect_true(r < 0.1)

  # Opposite directions cancel
  angles <- c(0, pi)
  expect_equal(circ_r(angles), 0, tolerance = 1e-10)
})

test_that("circ_r applies binning correction", {
  angles <- c(0, 0.1, -0.1)
  r_uncorrected <- circ_r(angles, d = 0)
  r_corrected <- circ_r(angles, d = 0.5)
  # Correction factor increases r

  expect_true(r_corrected > r_uncorrected)
})

test_that("circ_vtest detects non-uniformity with known direction", {
  set.seed(123)
  # Data concentrated around 0
  angles <- rnorm(50, mean = 0, sd = 0.3)
  result <- circ_vtest(angles, dir = 0)
  expect_true(result$pval < 0.05)
  expect_true(result$v > 0)

  # Data concentrated around pi, testing against dir = 0
  angles_pi <- rnorm(50, mean = pi, sd = 0.3)
  result_pi <- circ_vtest(angles_pi, dir = 0)
  # Should not be significant for wrong direction
  expect_true(result_pi$pval > result$pval)
})

test_that("circ_rayleigh detects non-uniformity", {
  set.seed(456)
  # Concentrated data
  concentrated <- rnorm(100, mean = 0, sd = 0.2)
  result <- circ_rayleigh(concentrated)
  expect_true(result$pval < 0.01)

  # Uniform data
  uniform <- seq(0, 2 * pi, length.out = 101)[1:100]
  result_uniform <- circ_rayleigh(uniform)
  expect_true(result_uniform$pval > 0.1)
})

test_that("fdr_bh correctly identifies significant tests", {
  pvals <- c(0.001, 0.01, 0.02, 0.03, 0.05, 0.1, 0.5)

  result <- fdr_bh(pvals, q = 0.05)
  expect_true(is.logical(result$h))
  expect_equal(length(result$h), length(pvals))
  expect_equal(length(result$adj_p), length(pvals))

  # First p-value should be significant
  expect_true(result$h[1])
  # Last p-value should not be significant
  expect_false(result$h[length(pvals)])

  # Adjusted p-values should be >= original
  expect_true(all(result$adj_p >= pvals))
})

test_that("fdr_bh handles edge cases",
{
  # All significant
  pvals_all_sig <- rep(0.001, 5)
  result <- fdr_bh(pvals_all_sig, q = 0.05)
  expect_true(all(result$h))

  # None significant
  pvals_none_sig <- rep(0.9, 5)
  result <- fdr_bh(pvals_none_sig, q = 0.05)
  expect_false(any(result$h))
  expect_equal(result$crit_p, 0)
})

# --- circ_mean matrix input ---

test_that("circ_mean works on matrix along dim=1", {
  mat <- matrix(c(0, pi/2, 0, pi/2), nrow = 2, ncol = 2)
  result <- circ_mean(mat, dim = 1)
  expect_length(result, 2)
  # Each column is (0, pi/2) -> mean should be pi/4
  expect_equal(as.numeric(result), rep(pi / 4, 2), tolerance = 1e-10)
})

test_that("circ_mean errors on weight dimension mismatch", {
  angles <- c(0, pi/2, pi)
  w <- c(1, 2)  # wrong length
  expect_error(circ_mean(angles, w = w), "dimensions do not match")
})

# --- circ_r matrix input ---

test_that("circ_r works on matrix with d correction", {
  mat <- matrix(c(0, 0.1, -0.1, 0, 0.1, -0.1), nrow = 3, ncol = 2)
  r <- circ_r(mat, d = 0.5, dim = 1)
  expect_length(r, 2)
  expect_true(all(r > 0))
})

# --- circ_vtest weight mismatch ---

test_that("circ_vtest errors on weight length mismatch", {
  angles <- c(0, pi/2, pi)
  w <- c(1, 2)
  expect_error(circ_vtest(angles, dir = 0, w = w), "dimensions do not match")
})

# --- fdr_bh single p-value ---

test_that("fdr_bh handles single p-value", {
  result <- fdr_bh(0.01, q = 0.05)
  expect_length(result$h, 1)
  expect_true(result$h[1])

  result2 <- fdr_bh(0.5, q = 0.05)
  expect_false(result2$h[1])
})
