test_that("oscillation_score increases with oscillation amplitude (synthetic simulation)", {
  set.seed(123)
  dt <- 0.002
  fs <- 1 / dt
  simulate_amp <- function(amp) {
    nsubj <- 5
    freq <- 6
    ddt <- 0.001
    Ttot <- 1
    out <- numeric(nsubj)
    for (s in seq_len(nsubj)) {
      timeline <- seq(0, Ttot, by = ddt)
      ratetrace <- (1 + amp * sin(2 * pi * freq * timeline)) * 15 / length(timeline)
      events <- stats::runif(length(timeline)) < ratetrace
      RTs <- timeline[events]
      mct <- make_continuous_trace(RTs, dt = dt, warn = FALSE)
      os <- oscillation_score(mct$signal, fs = fs, flim = c(1, 20), warnings = FALSE)
      sur <- oscillation_score_surrogates(mct$signal, fs = fs, flim = c(1, 20), nrep = 10, fpeak = os$fosc, warnings = FALSE)
      out[s] <- (log(os$oscore) - mean(log(sur$oscore_rp), na.rm = TRUE)) / sd(log(sur$oscore_rp), na.rm = TRUE)
    }
    out
  }
  z0 <- simulate_amp(0)
  z1 <- simulate_amp(0.8)
  expect_gt(mean(z1, na.rm = TRUE), mean(z0, na.rm = TRUE) + 1)
})
