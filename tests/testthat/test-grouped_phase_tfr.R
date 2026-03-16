make_group_burst <- function(n_subject = 6, n_bin = 120, fs = 30, freq = 6, noise_sd = 0.15) {
  t <- seq(0, by = 1 / fs, length.out = n_bin)
  env <- exp(-0.5 * ((seq_len(n_bin) - 60) / 10)^2)
  x <- vapply(seq_len(n_subject), function(i) {
    amp <- stats::runif(1, 0.9, 1.1)
    amp * env * sin(2 * pi * freq * t) + stats::rnorm(length(t), sd = noise_sd)
  }, numeric(n_bin))
  t(x)
}

make_phase_locked_pair <- function(n_subject = 6, n_bin = 100, fs = 25, freq = 5, phase_shift = pi / 2) {
  t <- seq(0, by = 1 / fs, length.out = n_bin)
  x1 <- vapply(seq_len(n_subject), function(i) {
    sin(2 * pi * freq * t) + 0.05 * stats::rnorm(length(t))
  }, numeric(n_bin))
  x2 <- vapply(seq_len(n_subject), function(i) {
    sin(2 * pi * freq * t + phase_shift) + 0.05 * stats::rnorm(length(t))
  }, numeric(n_bin))
  list(x1 = t(x1), x2 = t(x2))
}

test_that("group_tfr localizes burst power near target time and frequency", {
  set.seed(201)
  x <- make_group_burst()
  freqs <- c(4, 6, 8)

  res <- group_tfr(
    x,
    fs = 30,
    freqs = freqs,
    bandwidth = 1,
    reflect = TRUE,
    edge = 4
  )

  expect_equal(dim(res$observed$map), c(length(freqs), 120 - 8))
  peak_idx <- which(res$observed$map == max(res$observed$map), arr.ind = TRUE)[1, ]
  expect_equal(freqs[peak_idx[1]], 6)
  expect_true(abs(res$observed$time[peak_idx[2]] - 60) <= 12)
})

test_that("group_tfr_test returns finite statistics and clusters", {
  set.seed(202)
  x <- make_group_burst(n_subject = 7, n_bin = 110, fs = 30, freq = 6, noise_sd = 0.1)
  freqs <- c(4, 6, 8)

  res <- group_tfr_test(
    x,
    fs = 30,
    freqs = freqs,
    bandwidth = 1,
    null = "shuffle_labels",
    nrep = 25,
    reflect = TRUE,
    edge = 3,
    seed = 91
  )

  expect_equal(dim(res$null$map), c(25L, length(freqs), 110 - 6))
  expect_equal(dim(res$stats$z), c(length(freqs), 110 - 6))
  expect_true(all(is.finite(res$stats$p)))
  expect_true(is.list(res$stats$clusters))
  target_f <- which(freqs == 6)
  target_t <- which.min(abs(res$observed$time - 60))
  expect_true(is.finite(res$stats$z[target_f, target_t]))
})

test_that("group_phase_consistency detects strong phase alignment", {
  set.seed(203)
  pair <- make_phase_locked_pair()

  res <- group_phase_consistency(
    pair$x1,
    fs = 25,
    freq = 5,
    bandwidth = 0.75,
    reflect = TRUE,
    edge = 4
  )

  expect_equal(dim(res$observed$phase), c(6L, 100 - 8))
  expect_true(mean(res$stats$r, na.rm = TRUE) > 0.9)
  expect_true(all(res$stats$p >= 0 & res$stats$p <= 1))
})

test_that("group_phase_difference recovers known phase offset", {
  set.seed(204)
  pair <- make_phase_locked_pair(phase_shift = pi / 2)

  res <- group_phase_difference(
    pair$x1,
    pair$x2,
    fs = 25,
    freq = 5,
    bandwidth = 0.75,
    reflect = TRUE,
    edge = 4
  )

  mean_diff <- circ_mean(res$observed$mean_phase_difference)
  expect_true(abs(((mean_diff + pi / 2) + pi) %% (2 * pi) - pi) < 0.35)
  expect_true(mean(res$stats$r, na.rm = TRUE) > 0.85)
})

test_that("group_tfr and phase helpers validate arguments", {
  x <- matrix(rnorm(40), nrow = 4)

  expect_error(group_tfr(x, fs = 10, freqs = c(0, 2)), "strictly between 0 and Nyquist")
  expect_error(group_tfr(x, fs = 10, freqs = 2, edge = 5), "trims away all bins")
  expect_error(group_phase_difference(x, x[, -1, drop = FALSE], fs = 10, freq = 2), "matching bins")
})
