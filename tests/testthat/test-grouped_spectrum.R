make_group_sine <- function(n_subject = 6, n_bin = 120, fs = 30, freq = 6, noise_sd = 0.2) {
  t <- seq(0, by = 1 / fs, length.out = n_bin)
  x <- vapply(seq_len(n_subject), function(i) {
    phase <- stats::runif(1, -pi, pi)
    sin(2 * pi * freq * t + phase) + stats::rnorm(length(t), sd = noise_sd)
  }, numeric(n_bin))
  t(x)
}

test_that("group_spectrum recovers dominant group frequency", {
  set.seed(101)
  x <- make_group_sine()

  res <- group_spectrum(x, fs = 30, flim = c(3, 10), taper = "hann")

  expect_true(is.list(res))
  expect_true(all(c("observed", "meta") %in% names(res)))
  expect_equal(nrow(res$observed$per_subject), nrow(x))
  expect_equal(length(res$observed$freq), ncol(res$observed$per_subject))
  expect_true(abs(res$observed$peak - 6) <= 0.5)
})

test_that("group_spectrum accepts grouped-series objects", {
  set.seed(102)
  x <- make_group_sine(n_subject = 4, n_bin = 90, fs = 30, freq = 5)
  gx <- list(
    data = x,
    subjects = paste0("s", 1:4),
    bins = seq(200, by = 1000 / 30, length.out = 90),
    measure = "hr",
    meta = list(source = "test")
  )

  res <- group_spectrum(gx, fs = 30, flim = c(3, 8))

  expect_equal(res$meta$subjects, gx$subjects)
  expect_equal(res$meta$bins, gx$bins)
  expect_equal(res$meta$measure, gx$measure)
  expect_equal(res$meta$source, "test")
})

test_that("group_spectrum_test detects oscillatory structure under shuffle null", {
  set.seed(103)
  x <- make_group_sine(n_subject = 8, n_bin = 150, fs = 30, freq = 6, noise_sd = 0.1)

  res <- group_spectrum_test(
    x,
    fs = 30,
    flim = c(3, 10),
    null = "shuffle_labels",
    nrep = 80,
    taper = "hann",
    p_adjust = "fdr",
    seed = 55
  )

  idx <- which.min(abs(res$observed$freq - 6))
  expect_equal(dim(res$null$power), c(80L, length(res$observed$freq)))
  expect_true(is.finite(res$stats$z[idx]))
  expect_lt(res$stats$p_adj[idx], 0.1)
})

test_that("group_spectrum_test circular shift returns finite outputs", {
  set.seed(104)
  x <- make_group_sine(n_subject = 5, n_bin = 100, fs = 25, freq = 4, noise_sd = 0.3)

  res <- group_spectrum_test(
    x,
    fs = 25,
    flim = c(2, 8),
    null = "circular_shift",
    nrep = 20,
    seed = 77
  )

  expect_equal(length(res$stats$p), length(res$observed$freq))
  expect_equal(length(res$stats$significant), length(res$observed$freq))
  expect_true(all(is.finite(res$null$power)))
})

test_that("group_spectrum validates inputs", {
  expect_error(group_spectrum(matrix(1, nrow = 1, ncol = 1), fs = 10), "at least one subject and two bins")
  expect_error(group_spectrum(matrix(c(1, NA), nrow = 1), fs = 10), "finite")
  expect_error(group_spectrum(matrix(1:6, nrow = 2), fs = -1), "positive scalar")
  expect_error(group_spectrum_test(matrix(1:6, nrow = 2), fs = 10, nrep = 0), "positive scalar")
})
