test_that("pairwise_phase_consistency computes expected PPC", {
  angles <- matrix(c(0, 0, pi / 2, pi / 2), nrow = 2, byrow = TRUE)
  ppc <- pairwise_phase_consistency(angles, dim = 1)
  expect_true(all(ppc >= 0))
})

test_that("pairwise_phase_consistency is 1 for identical phases", {
  angles <- matrix(0, nrow = 3, ncol = 5)
  ppc <- pairwise_phase_consistency(angles, dim = 1)
  expect_true(all(abs(ppc - 1) < 1e-8))
})

test_that("pairwise_phase_consistency ignores NA values like MATLAB isnan", {
  angles <- matrix(c(0, 0,
                     0, NA,
                     0, 0),
                   nrow = 3, byrow = TRUE)
  ppc <- pairwise_phase_consistency(angles, dim = 1)
  expect_true(all(is.finite(ppc)))
  expect_true(all(abs(ppc - 1) < 1e-8))
})

test_that("u_score_matrix returns Z and significance", {
  dat <- c(1, 2, 3)
  ref <- rbind(c(0, 0, 0), c(0, 0, 0))
  res <- u_score_matrix(dat, ref, alpha = 0.05)
  expect_equal(length(res$Z), length(dat))
  expect_true(all(res$sgnf == 1))
})

test_that("paired_tscore handles one-sample", {
  dat <- matrix(c(1, 2, 3, 4), nrow = 2)
  res <- paired_tscore(dat, dim = 1)
  expect_true(all(dim(res$t) == c(1, 2)))
})

test_that("nonparam_pval clamps negative p", {
  data <- c(5)
  ref <- c(1, 1, 1)
  res <- nonparam_pval(data, ref, alpha = 0.05)
  expect_true(res$p >= 0)
  expect_true(res$h == 1)
})

# --- paired_tscore with dat2 (paired) ---

test_that("paired_tscore computes paired difference", {
  dat1 <- matrix(c(5, 6, 7, 8), nrow = 2)
  dat2 <- matrix(c(1, 2, 3, 4), nrow = 2)
  res <- paired_tscore(dat1, dat2, dim = 1)
  expect_true(all(dim(res$t) == c(1, 2)))
  # Differences are all positive, so t should be positive
  expect_true(all(res$t > 0, na.rm = TRUE))
})

test_that("paired_tscore returns NA for zero SD", {
  # All values equal along dim → SD = 0
  dat <- matrix(c(5, 5, 3, 3), nrow = 2)
  res <- paired_tscore(dat, dim = 1)
  expect_true(all(is.na(res$t)))
})

# --- nonparam_pval additional branches ---

test_that("nonparam_pval with value below median", {
  ref <- c(5, 6, 7, 8, 9, 10)
  res <- nonparam_pval(c(4), ref, alpha = 0.05)
  expect_true(res$p >= 0)
  expect_true(res$p <= 1)
})

test_that("nonparam_pval with value above median", {
  ref <- c(1, 2, 3, 4, 5)
  res <- nonparam_pval(c(10), ref, alpha = 0.05)
  expect_true(res$p >= 0)
  expect_true(res$h == 1)
})
