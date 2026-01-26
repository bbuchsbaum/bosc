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

# --- oscillation_score additional branches ---

test_that("oscillation_score warns on sum < 3", {
  sig <- c(1, 0, 0, 0, 0)
  expect_warning(
    res <- oscillation_score(sig, fs = 100, flim = c(1, 10), warnings = TRUE),
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
  res <- oscillation_score(sig, fs = fs, flim = c(1, 50),
                           fcor = TRUE, warnings = FALSE)
  expect_true(is.numeric(res$oscore))
})

test_that("oscillation_score validates inputs", {
  expect_error(oscillation_score("abc", fs = 100, flim = c(1, 10)), "numeric")
  expect_error(oscillation_score(1:10, fs = -1, flim = c(1, 10)), "positive")
  expect_error(oscillation_score(1:10, fs = 100, flim = c(1)), "length 2")
})

# --- oscillation_score_surrogates validation ---

test_that("oscillation_score_surrogates validates inputs", {
  expect_error(oscillation_score_surrogates(1:10, fs = 100, flim = c(1, 10),
                                             nrep = 0, fpeak = 5), "positive")
  expect_error(oscillation_score_surrogates(1:10, fs = -1, flim = c(1, 10),
                                             nrep = 5, fpeak = 5), "positive")
  expect_error(oscillation_score_surrogates(1:10, fs = 100, flim = c(1, 10),
                                             nrep = 5, fpeak = NULL), "fpeak must be provided")
})
