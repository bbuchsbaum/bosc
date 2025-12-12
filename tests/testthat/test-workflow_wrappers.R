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
  res <- oscillation_score_z(sig, fs = fs, flim = c(1, 20), nrep = 25)
  valid <- res$surrogates$oscore_rp
  valid <- valid[is.finite(valid) & valid > 0]
  z_manual <- (log(res$oscore) - mean(log(valid))) / sd(log(valid))
  expect_equal(res$z, z_manual)
})

test_that("oscillation_score_z tidy=TRUE returns a data.frame", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  tidy_res <- oscillation_score_z(sig, fs = fs, flim = c(1, 20), nrep = 10, tidy = TRUE)
  expect_s3_class(tidy_res, "data.frame")
  expect_true(all(c("oscore", "fosc", "z", "pval", "significant") %in% names(tidy_res)))
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
