test_that("make_continuous_trace maps events and keeps length", {
  ev <- c(0, 0.1, 0.2)
  res <- make_continuous_trace(ev, dt = 0.1, warn = FALSE)
  expect_equal(res$tspan, c(0, 0.1, 0.2, 0.3))
  expect_equal(res$signal, c(1, 1, 1, 0))
})

test_that("make_continuous_trace trims quantiles", {
  ev <- c(0, 1, 2, 3)
  res <- make_continuous_trace(ev, dt = 1, quantlim = c(0.25, 0.75), warn = FALSE)
  expect_equal(res$tspan, c(0, 1, 2))
  expect_equal(res$signal, c(1, 1, 1))
})

test_that("autocorr_centered matches manual convolution", {
  x <- c(1, 2)
  expect_equal(autocorr_centered(x), c(2, 5, 2))
})

test_that("spectral_peak identifies dominant frequency", {
  fs <- 100
  t <- seq(0, 1, by = 1 / fs)
  data <- sin(2 * pi * 10 * t)
  sp <- spectral_peak(data, fs = fs, flim = c(1, 20))
  expect_true(abs(sp$freq - 10) < 0.5)
  expect_length(sp$fxx, length(sp$spectrum))
})
