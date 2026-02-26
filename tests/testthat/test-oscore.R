test_that("oscillation_score detects a simple periodic spike train", {
  fs <- 1000
  signal <- numeric(fs)
  signal[seq(100, 900, 100)] <- 1
  res <- oscillation_score(signal, fs = fs, flim = c(1, 50), warnings = FALSE)
  expect_true(is.finite(res$oscore))
  expect_true(abs(res$fosc - 10) < 1.5) # allow some tolerance
})

test_that("oscillation_score returns NA for empty-like input", {
  res <- oscillation_score(rep(0, 100), fs = 100, flim = c(1, 10), warnings = FALSE)
  expect_true(is.na(res$oscore))
})

test_that("oscillation_score_surrogates returns replicates", {
  fs <- 100
  signal <- numeric(500)
  signal[seq(50, 450, 50)] <- 1
  sur <- oscillation_score_surrogates(signal,
                                      fs = fs,
                                      flim = c(1, 20),
                                      nrep = 5,
                                      fpeak = 2,
                                      warnings = FALSE)
  expect_length(sur$oscore_rp, 5)
  expect_true(any(is.finite(sur$oscore_rp)))
})

test_that("oscillation_score_surrogates supports phase-randomized continuous surrogates", {
  set.seed(5)
  fs <- 200
  t <- seq(0, 2, by = 1 / fs)
  signal <- sin(2 * pi * 8 * t) + 0.2 * rnorm(length(t))
  sur <- oscillation_score_surrogates(
    signal,
    fs = fs,
    flim = c(4, 20),
    nrep = 4,
    fpeak = 8,
    surrogate_method = "phase_randomized",
    signal_mode = "continuous",
    warnings = FALSE
  )
  expect_equal(sur$surrogate_method, "phase_randomized")
  expect_equal(sur$signal_mode, "continuous")
  expect_length(sur$oscore_rp, 4)
  expect_true(all(vapply(sur$signrep, length, integer(1)) == length(signal)))
  expect_true(any(is.finite(sur$oscore_rp)))

  # Phase randomization should preserve the FFT amplitude spectrum.
  amp0 <- Mod(fft(signal - mean(signal)))
  amp1 <- Mod(fft(sur$signrep[[1]] - mean(sur$signrep[[1]])))
  expect_equal(amp1, amp0, tolerance = 1e-8)
})

# --- oscillation_score additional branches ---

test_that("oscillation_score warns on sum < 3", {
  sig <- c(1, 0, 0, 0, 0)
  expect_warning(
    res <- oscillation_score(sig, fs = 100, flim = c(1, 10), warnings = TRUE, signal_mode = "event"),
    "nearly.*empty"
  )
  expect_true(is.na(res$oscore))
})

test_that("oscillation_score uses quantlim when specified", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  res <- oscillation_score(sig, fs = fs, flim = c(1, 20),
                           quantlim = c(0.1, 0.9), warnings = FALSE)
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score works with smoothach=FALSE", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  res <- oscillation_score(sig, fs = fs, flim = c(1, 50),
                           smoothach = FALSE, warnings = FALSE)
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score uses fpeak override", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  res <- oscillation_score(sig, fs = fs, flim = c(1, 50),
                           fpeak = 15, warnings = FALSE)
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score accepts taper='hann'", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  res <- oscillation_score(sig, fs = fs, flim = c(1, 50),
                           taper = "hann", warnings = FALSE)
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score accepts fcor=TRUE", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  expect_silent(
    res <- oscillation_score(sig, fs = fs, flim = c(1, 50),
                             fcor = TRUE, warnings = FALSE)
  )
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score validates inputs", {
  expect_error(oscillation_score("abc", fs = 100, flim = c(1, 10)), "numeric")
  expect_error(oscillation_score(1:10, fs = -1, flim = c(1, 10)), "positive")
  expect_error(oscillation_score(1:10, fs = 100, flim = c(1)), "length 2")
  expect_error(oscillation_score(1:10, fs = 100, flim = c(NA, 10)), "finite")
  expect_error(oscillation_score(1:10, fs = 100, flim = c(0, 10)), "increasing positive")
})

test_that("oscillation_score handles continuous oscillatory input in auto mode", {
  set.seed(11)
  fs <- 1000
  t <- seq(0, 2, by = 1 / fs)
  sig <- sin(2 * pi * 10 * t) + rnorm(length(t), sd = 0.2)
  res <- oscillation_score(sig, fs = fs, flim = c(5, 20), warnings = FALSE)
  expect_true(is.finite(res$oscore))
  expect_true(is.finite(res$fosc))
  expect_true(res$fosc >= 5 && res$fosc <= 20)
})

test_that("oscillation_score reports near-zero variance for continuous mode", {
  sig <- rep(0.5, 500)
  expect_warning(
    res <- oscillation_score(sig, fs = 100, flim = c(5, 20), signal_mode = "continuous", warnings = TRUE),
    "near-zero variance"
  )
  expect_true(is.na(res$oscore))
})

test_that("oscillation_score warns when effective band is invalid", {
  fs <- 100
  sig <- numeric(100)
  sig[c(10, 50, 90)] <- 1
  expect_warning(
    res <- oscillation_score(sig, fs = fs, flim = c(5, 40), mincycles = 80, warnings = TRUE),
    "Effective frequency bounds are invalid"
  )
  expect_true(is.na(res$oscore))
  expect_length(res$flim, 2)
})

test_that("oscillation_score enforces minfreqbandwidth", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 10)] <- 1
  expect_warning(
    res <- oscillation_score(
      sig, fs = fs, flim = c(3, 20), minfreqbandwidth = 25, warnings = TRUE
    ),
    "minimal frequency bandwidth"
  )
  expect_true(is.na(res$oscore))
})

test_that("oscillation_score allows zero peakwind and optional plotting", {
  skip_if_not_installed("ggplot2")
  fs <- 500
  sig <- numeric(fs)
  sig[seq(50, 450, 50)] <- 1
  expect_silent(
    res <- oscillation_score(
      sig, fs = fs, flim = c(3, 20), peakwind = 0, plot = TRUE, warnings = FALSE
    )
  )
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score returns NA when peak frequency is unavailable", {
  call_n <- 0L
  local_mocked_bindings(
    spectral_peak = function(...) {
      call_n <<- call_n + 1L
      if (call_n == 1L) {
        list(freq = NA_real_, fxx = c(5, 6), spectrum = c(1, 2))
      } else {
        list(freq = NA_real_, fxx = c(5, 6), spectrum = c(1, 2))
      }
    },
    .package = "bosc"
  )
  sig <- c(0, 3, 0)
  res <- oscillation_score(sig, fs = 100, flim = c(1, 10), mincycles = 0, signal_mode = "event", warnings = FALSE)
  expect_true(is.na(res$oscore))
})

test_that("oscillation_score handles non-finite peak index and zero-mean spectrum", {
  # First branch: index resolution fails (freq axis all NA).
  call_n <- 0L
  local_mocked_bindings(
    spectral_peak = function(...) {
      call_n <<- call_n + 1L
      if (call_n == 1L) {
        list(freq = 5, fxx = c(NA_real_), spectrum = c(1))
      } else {
        list(freq = 5, fxx = c(1), spectrum = c(1))
      }
    },
    .package = "bosc"
  )
  sig <- c(0, 3, 0)
  res_idx <- suppressWarnings(
    oscillation_score(sig, fs = 100, flim = c(1, 10), mincycles = 0, signal_mode = "event", warnings = FALSE)
  )
  expect_true(is.na(res_idx$oscore))

  # Second branch: denominator mean is invalid.
  call_n <- 0L
  local_mocked_bindings(
    spectral_peak = function(...) {
      call_n <<- call_n + 1L
      if (call_n == 1L) {
        list(freq = 5, fxx = c(4, 5, 6), spectrum = c(1, 2, 1))
      } else {
        list(freq = 5, fxx = c(4, 5, 6), spectrum = c(0, 0, 0))
      }
    },
    .package = "bosc"
  )
  expect_warning(
    res_den <- oscillation_score(sig, fs = 100, flim = c(1, 10), mincycles = 0, signal_mode = "event", warnings = TRUE),
    "Mean spectrum is zero/invalid"
  )
  expect_true(is.na(res_den$oscore))
})

# --- oscillation_score_surrogates validation ---

test_that("oscillation_score_surrogates validates inputs", {
  expect_error(oscillation_score_surrogates(1:10, fs = 100, flim = c(1, 10),
                                             nrep = 0, fpeak = 5), "positive")
  expect_error(oscillation_score_surrogates(1:10, fs = -1, flim = c(1, 10),
                                             nrep = 5, fpeak = 5), "positive")
  expect_error(oscillation_score_surrogates(1:10, fs = 100, flim = c(1, 10),
                                             nrep = 5, fpeak = NULL), "fpeak must be provided")
  expect_error(oscillation_score_surrogates(1:10, fs = 100, flim = c(NA, 10),
                                             nrep = 5, fpeak = 5), "increasing positive")
  expect_error(oscillation_score_surrogates(1:10, fs = 100, flim = c(10, 1),
                                             nrep = 5, fpeak = 5), "increasing positive")
})

test_that("oscillation_score_surrogates warns when keep_trend is ignored", {
  fs <- 100
  sig <- numeric(300)
  sig[seq(30, 270, 30)] <- 1
  expect_warning(
    res <- oscillation_score_surrogates(
      sig,
      fs = fs,
      flim = c(2, 20),
      nrep = 3,
      fpeak = 5,
      keep_trend = TRUE,
      surrogate_method = "phase_randomized",
      warnings = TRUE
    ),
    "keep_trend ignored"
  )
  expect_equal(res$surrogate_method, "phase_randomized")
})

test_that("oscillation_score_surrogates early returns for sparse or flat signals", {
  expect_warning(
    sparse <- oscillation_score_surrogates(
      c(1, 0, 0, 0),
      fs = 100,
      flim = c(1, 10),
      nrep = 2,
      fpeak = 4,
      signal_mode = "event",
      warnings = TRUE
    ),
    "could not be computed"
  )
  expect_true(is.na(sparse$oscore_rp[1]))

  expect_warning(
    flat <- oscillation_score_surrogates(
      rep(0.5, 50),
      fs = 100,
      flim = c(1, 10),
      nrep = 2,
      fpeak = 4,
      signal_mode = "continuous",
      warnings = TRUE
    ),
    "near-zero variance"
  )
  expect_true(is.na(flat$oscore_rp[1]))
})

test_that("oscillation_score_surrogates trend mode supports fallback and fitted path", {
  set.seed(101)
  fs <- 100
  sig <- numeric(600)
  sig[sample(50:550, 140)] <- 1

  expect_warning(
    fallback <- oscillation_score_surrogates(
      sig,
      fs = fs,
      flim = c(1, 20),
      nrep = 2,
      fpeak = 6,
      keep_trend = TRUE,
      surrogate_method = "auto",
      trend_dist = "notadist",
      warnings = TRUE
    ),
    "Cannot fit distribution"
  )
  expect_equal(fallback$surrogate_method, "event_shuffle")

  fitted <- oscillation_score_surrogates(
    sig,
    fs = fs,
    flim = c(1, 20),
    nrep = 3,
    fpeak = 6,
    surrogate_method = "event_trend",
    trend_dist = "gamma",
    trend_alpha = -1,
    warnings = FALSE
  )
  expect_equal(fitted$surrogate_method, "event_trend")
  expect_false(is.null(fitted$trendfit$pd))
  expect_true(any(vapply(fitted$signrep, sum, numeric(1)) > 0))
})

test_that("phase_randomize_signal handles short and centered-zero inputs", {
  phase_rand <- getFromNamespace("phase_randomize_signal", "bosc")
  expect_equal(phase_rand(c(1, 2)), c(1, 2))
  expect_equal(phase_rand(rep(3, 8)), rep(3, 8))
})
