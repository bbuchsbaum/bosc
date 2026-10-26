simulate_events <- function(amp, freq = 6, dur = 10, rate = 100, ddt = 0.001) {
  timeline <- seq(0, dur, by = ddt)
  ratetrace <- (1 + amp * sin(2 * pi * freq * timeline)) * rate * ddt
  timeline[stats::runif(length(timeline)) < ratetrace]
}

test_that("oscillation_score recovers a rhythmic event train (synthetic simulation)", {
  set.seed(123)
  fs <- 500
  score_for <- function(amp) {
    vapply(seq_len(4), function(i) {
      mct <- make_continuous_trace(simulate_events(amp), dt = 1 / fs, warn = FALSE)
      os <- oscillation_score(mct$signal, fs = fs, flim = c(1, 20), warnings = FALSE)
      c(oscore = os$oscore, fosc = os$fosc)
    }, numeric(2))
  }
  rhythmic <- score_for(0.8)
  flat <- score_for(0)

  expect_true(all(abs(rhythmic["fosc", ] - 6) < 1))
  expect_gt(min(rhythmic["oscore", ]), 3 * max(flat["oscore", ]))
})

test_that("surrogate z-scores separate rhythmic from flat event trains", {
  set.seed(42)
  fs <- 500
  z_for <- function(amp) {
    mct <- make_continuous_trace(simulate_events(amp), dt = 1 / fs, warn = FALSE)
    os <- oscillation_score(mct$signal, fs = fs, flim = c(1, 20), warnings = FALSE)
    sur <- oscillation_score_surrogates(
      mct$signal, fs = fs, flim = c(1, 20), nrep = 20, fpeak = os$fosc, warnings = FALSE
    )
    (log(os$oscore) - mean(log(sur$oscore_rp), na.rm = TRUE)) /
      stats::sd(log(sur$oscore_rp), na.rm = TRUE)
  }
  expect_gt(z_for(0.8), 8)
  expect_lt(z_for(0), 6)
})

test_that("make_continuous_trace output is exactly non-negative and detected as event-like", {
  set.seed(1)
  mct <- make_continuous_trace(sort(stats::runif(200, 0, 5)), dt = 0.002, warn = FALSE)
  expect_true(all(mct$signal >= 0))
  sur <- oscillation_score_surrogates(
    mct$signal, fs = 500, flim = c(1, 20), nrep = 2, fpeak = 6, warnings = FALSE
  )
  expect_equal(sur$signal_mode, "event")
  expect_equal(sur$surrogate_method, "event_shuffle")
})
