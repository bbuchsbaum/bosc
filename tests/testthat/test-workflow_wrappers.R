test_that("oscillation_score_config matches direct call", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  config <- list(fs = fs, flim = c(1, 50), warnings = "off", taper = "hanning")
  res_cfg <- oscillation_score_config(config, sig)
  res_dir <- oscillation_score(sig, fs = fs, flim = c(1, 50), warnings = FALSE, taper = "hanning")
  expect_equal(res_cfg$oscore, res_dir$oscore)
  expect_equal(res_cfg$fosc, res_dir$fosc)
})

test_that("oscillation_score_z computes log-Z using returned surrogates", {
  set.seed(1)
  fs <- 500
  sig <- numeric(fs)
  sig[seq(50, 450, 50)] <- 1
  res <- oscillation_score_z(sig, fs = fs, flim = c(1, 20), nrep = 25, ci_nboot = 0)
  valid <- res$surrogates$oscore_rp
  valid <- valid[is.finite(valid) & valid > 0]
  z_manual <- (log(res$oscore) - mean(log(valid))) / sd(log(valid))
  expect_equal(res$z, z_manual)
})

test_that("oscillation_score_z tidy=TRUE returns a data.frame", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  tidy_res <- oscillation_score_z(sig, fs = fs, flim = c(1, 20), nrep = 10, ci_nboot = 0, tidy = TRUE)
  expect_s3_class(tidy_res, "data.frame")
  expect_true(all(c("oscore", "fosc", "z", "pval", "significant") %in% names(tidy_res)))
})

test_that("oscillation_score_z returns bootstrap CI fields", {
  set.seed(9)
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  res <- oscillation_score_z(sig, fs = fs, flim = c(1, 20), nrep = 20, ci_nboot = 100)
  expect_length(res$z_ci, 2)
  expect_equal(res$ci_level, 0.95)
  expect_equal(res$ci_nboot, 100)
})

test_that("oscillation_score_z supports phase-randomized continuous mode", {
  set.seed(12)
  fs <- 400
  t <- seq(0, 2, by = 1 / fs)
  sig <- sin(2 * pi * 9 * t) + rnorm(length(t), sd = 0.2)
  res <- oscillation_score_z(
    sig,
    fs = fs,
    flim = c(5, 20),
    nrep = 15,
    surrogate_method = "phase_randomized",
    signal_mode = "continuous",
    ci_nboot = 50
  )
  expect_true(is.finite(res$oscore))
  expect_true(is.finite(res$fosc))
  expect_equal(res$surrogates$surrogate_method, "phase_randomized")
})

test_that("phase_at_events returns phases aligned with events", {
  events <- seq(0.05, 0.95, by = 0.1)  # ~10 Hz
  ph <- phase_at_events(events, dt = 0.001, fosc = 10, bandwidth = 1)
  expect_length(ph, length(events))
  expect_true(all(is.finite(ph)))

  events2 <- c(events, NA_real_)
  ph2 <- phase_at_events(events2, dt = 0.001, fosc = 10, bandwidth = 1)
  expect_true(is.na(ph2[length(ph2)]))
})

# --- oscillation_score_config additional branches ---

test_that("oscillation_score_config errors on NULL config", {
  expect_error(oscillation_score_config(NULL, 1:10), "config must be a list")
})

test_that("oscillation_score_config converts character warnings", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  config <- list(fs = fs, flim = c(1, 50), warnings = "off", taper = "hanning")
  res <- oscillation_score_config(config, sig)
  expect_true(is.numeric(res$oscore))
})

# --- oscillation_score_surrogates_config ---

test_that("oscillation_score_surrogates_config forwards correctly", {
  fs <- 100
  sig <- numeric(500)
  sig[seq(50, 450, 50)] <- 1
  config <- list(fs = fs, flim = c(1, 20), nrep = 3, fpeak = 2, warnings = "off")
  res <- oscillation_score_surrogates_config(config, sig)
  expect_length(res$oscore_rp, 3)
})

test_that("oscillation_score_surrogates_config forwards surrogate method args", {
  set.seed(33)
  fs <- 200
  t <- seq(0, 1, by = 1 / fs)
  sig <- sin(2 * pi * 8 * t) + 0.1 * rnorm(length(t))
  config <- list(
    fs = fs,
    flim = c(4, 20),
    nrep = 3,
    fpeak = 8,
    surrogate_method = "phase_randomized",
    signal_mode = "continuous",
    warnings = "off"
  )
  res <- oscillation_score_surrogates_config(config, sig)
  expect_equal(res$surrogate_method, "phase_randomized")
  expect_equal(res$signal_mode, "continuous")
})

test_that("oscillation_score_surrogates_config errors on NULL config", {
  expect_error(oscillation_score_surrogates_config(NULL, 1:10), "config must be a list")
})

# --- deprecated wrappers ---

test_that("oscillation_score_cfg emits deprecation warning", {
  fs <- 1000
  sig <- numeric(fs)
  sig[seq(100, 900, 100)] <- 1
  cfg <- list(fs = fs, flim = c(1, 50), warnings = "off")
  expect_warning(oscillation_score_cfg(cfg, sig), "deprecated", ignore.case = TRUE)
})

test_that("oscillation_score_stats emits deprecation warning", {
  fs <- 100
  sig <- numeric(500)
  sig[seq(50, 450, 50)] <- 1
  cfg <- list(fs = fs, flim = c(1, 20), nrep = 3, fpeak = 2, warnings = "off")
  expect_warning(oscillation_score_stats(cfg, sig), "deprecated", ignore.case = TRUE)
})

test_that("oscore_z emits deprecation warning", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  expect_warning(
    oscore_z(sig, fs = fs, flim = c(1, 20), nrep = 5),
    "deprecated", ignore.case = TRUE
  )
})

# --- oscillation_score_z early return for invalid oscore ---

test_that("oscillation_score_z returns NA z for all-zero signal", {
  res <- oscillation_score_z(rep(0, 100), fs = 100, flim = c(1, 10), nrep = 5, ci_nboot = 0)
  expect_true(is.na(res$z))
  expect_true(is.na(res$pval))
  expect_false(res$significant)
  expect_null(res$surrogates)
})

test_that("oscillation_score_z validates CI arguments", {
  sig <- c(0, 1, 0, 1, 0, 1)
  expect_error(
    oscillation_score_z(sig, fs = 10, flim = c(1, 4), nrep = 5, ci_nboot = -1),
    "non-negative"
  )
  expect_error(
    oscillation_score_z(sig, fs = 10, flim = c(1, 4), nrep = 5, ci_level = 1.2),
    "in \\(0, 1\\)"
  )
})

# --- phase_at_events additional branches ---

test_that("phase_at_events with leave_one_out=TRUE", {
  events <- seq(0.05, 0.95, by = 0.1)
  ph <- phase_at_events(events, dt = 0.001, fosc = 10, bandwidth = 1,
                         leave_one_out = TRUE)
  expect_length(ph, length(events))
})

test_that("phase_at_events returns all NA for all-NA events", {
  events <- c(NA_real_, NA_real_, NA_real_)
  ph <- phase_at_events(events, dt = 0.001, fosc = 10, bandwidth = 1)
  expect_true(all(is.na(ph)))
})

test_that("phase_at_events validates inputs", {
  expect_error(phase_at_events("bad", dt = 0.001, fosc = 10), "events must be numeric")
  expect_error(phase_at_events(c(0.1, 0.2), dt = -1, fosc = 10), "positive scalar")
})
