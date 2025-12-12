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
