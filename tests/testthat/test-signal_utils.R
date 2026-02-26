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

test_that("autocorr_centered fft matches direct across sizes", {
  set.seed(99)
  for (n in c(5, 11, 257, 300)) {
    x <- rnorm(n)
    expect_equal(
      autocorr_centered(x, method = "fft"),
      autocorr_centered(x, method = "direct"),
      tolerance = 1e-8
    )
  }
})

test_that("autocorr_centered auto path tracks size heuristic", {
  set.seed(100)
  x_short <- rnorm(64)
  x_long <- rnorm(400)
  expect_equal(
    autocorr_centered(x_short, method = "auto"),
    autocorr_centered(x_short, method = "direct"),
    tolerance = 1e-10
  )
  expect_equal(
    autocorr_centered(x_long, method = "auto"),
    autocorr_centered(x_long, method = "fft"),
    tolerance = 1e-10
  )
})

test_that("spectral_peak identifies dominant frequency", {
  fs <- 100
  t <- seq(0, 1, by = 1 / fs)
  data <- sin(2 * pi * 10 * t)
  sp <- spectral_peak(data, fs = fs, flim = c(1, 20))
  expect_true(abs(sp$freq - 10) < 0.5)
  expect_length(sp$fxx, length(sp$spectrum))
})

# --- make_continuous_trace additional branches ---

test_that("make_continuous_trace drops events at remove_val with warning", {
  ev <- c(0, 0.1, 0.5, 0.2)
  expect_warning(
    res <- make_continuous_trace(ev, dt = 0.1, remove_val = 0.5, warn = TRUE),
    "Dropped"
  )
  # The 0.5 event should be removed
  expect_true(length(res$signal) > 0)
})

test_that("make_continuous_trace applies sd_smooth", {
  ev <- c(0, 0.1, 0.2)
  res <- make_continuous_trace(ev, dt = 0.01, sd_smooth = 0.05, warn = FALSE)
  # Smoothed signal should not be pure 0/1
  expect_true(any(res$signal > 0 & res$signal < 1))
})

test_that("make_continuous_trace applies width_block", {
  ev <- c(0, 0.1, 0.2)
  res <- make_continuous_trace(ev, dt = 0.01, width_block = 0.05, warn = FALSE)
  expect_true(length(res$signal) > 0)
})

test_that("make_continuous_trace returns empty for no events", {
  res <- make_continuous_trace(numeric(0), dt = 0.1, warn = FALSE)
  expect_length(res$signal, 0)
  expect_length(res$tspan, 0)
})

# --- spectral_peak additional branches ---

test_that("spectral_peak uses hann taper", {
  fs <- 100
  t <- seq(0, 1, by = 1 / fs)
  data <- sin(2 * pi * 10 * t)
  sp <- spectral_peak(data, fs = fs, flim = c(1, 20), taper = "hann")
  expect_true(is.finite(sp$freq))
})

test_that("spectral_peak applies fcor=TRUE", {
  fs <- 100
  t <- seq(0, 1, by = 1 / fs)
  data <- sin(2 * pi * 10 * t) + rnorm(length(t), sd = 0.1)
  expect_silent(sp <- spectral_peak(data, fs = fs, flim = c(1, 20), fcor = TRUE))
  expect_true(is.numeric(sp$freq))
})

test_that("spectral_peak fcor handles full band including 0 Hz", {
  fs <- 100
  t <- seq(0, 1, by = 1 / fs)
  data <- sin(2 * pi * 10 * t) + 0.25
  expect_silent(sp <- spectral_peak(data, fs = fs, flim = NULL, fcor = TRUE))
  expect_equal(length(sp$fxx), length(sp$spectrum))
  expect_true(any(is.finite(sp$spectrum)))
})

test_that("spectral_peak with flim=NULL returns full spectrum", {
  fs <- 100
  t <- seq(0, 1, by = 1 / fs)
  data <- sin(2 * pi * 10 * t)
  sp <- spectral_peak(data, fs = fs, flim = NULL)
  expect_true(length(sp$fxx) > 0)
})

test_that("spectral_peak validates fs", {
  expect_error(spectral_peak(1:10, fs = -1, flim = c(1, 10)), "positive scalar")
})
