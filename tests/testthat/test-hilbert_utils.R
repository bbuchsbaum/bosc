test_that("narrowband_hilbert recovers phase of a pure tone", {
  fs <- 200
  t <- seq(0, 1, by = 1 / fs)
  freq <- 10
  sig <- sin(2 * pi * freq * t)
  res <- narrowband_hilbert(sig, fs = fs, freqlim = c(8, 12), tol = 1e3)
  amp <- Mod(res$analytic)
  expect_true(mean(amp, na.rm = TRUE) > 0.5)
  sp <- spectral_peak(res$filtered, fs = fs, flim = c(1, 30))
  expect_true(abs(sp$freq - freq) < 1)
})

test_that("narrowband_hilbert flags instability", {
  fs <- 100
  sig <- c(rep(1, 10), rep(1000, 10), rep(1, 10))
  expect_warning(
    res <- narrowband_hilbert(sig, fs = fs, freqlim = c(0, 40), tol = 0.5),
    "Filter is unstable",
    fixed = TRUE
  )
  expect_true(any(is.na(res$analytic)))
})

# --- narrowband_hilbert additional branches ---

test_that("narrowband_hilbert uses lowpass when freqlim[1] <= 0", {
  fs <- 200
  t <- seq(0, 1, by = 1 / fs)
  sig <- sin(2 * pi * 5 * t)
  res <- narrowband_hilbert(sig, fs = fs, freqlim = c(0, 10), tol = 1e3)
  expect_true(length(res$filtered) == length(sig))
  expect_true(any(is.finite(res$analytic)))
})

test_that("narrowband_hilbert with demean=FALSE", {
  fs <- 200
  t <- seq(0, 1, by = 1 / fs)
  sig <- sin(2 * pi * 10 * t) + 5  # offset signal
  res <- narrowband_hilbert(sig, fs = fs, freqlim = c(8, 12), tol = 1e3, demean = FALSE)
  expect_true(length(res$filtered) == length(sig))
})

test_that("narrowband_hilbert validates inputs", {
  expect_error(narrowband_hilbert(1:10, fs = -1, freqlim = c(1, 10)), "positive scalar")
  expect_error(narrowband_hilbert(1:10, fs = 100, freqlim = c(1)), "length 2")
  expect_error(narrowband_hilbert(1:10, fs = 100, freqlim = c(1, 10), filtorder = -1), "positive")
})
