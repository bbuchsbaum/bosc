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
