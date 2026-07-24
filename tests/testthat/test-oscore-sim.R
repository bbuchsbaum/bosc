test_that("oscillation_score increases with oscillation amplitude (synthetic simulation)", {
  # This simulation validates the core contract of oscillation_score: a spike
  # train whose rate is modulated at a known frequency should score higher, and
  # recover that frequency, relative to an unmodulated (constant-rate) train.
  #
  # The check is deliberately over-powered rather than marginal: enough subjects
  # to shrink the between-subject SE, a long enough window, and a healthy event
  # rate. A previous version compared a surrogate-normalised z with only 5
  # subjects / 15 events and failed deterministically because the effect,
  # though real, was swamped by sampling noise.
  set.seed(123)
  dt <- 0.002
  fs <- 1 / dt
  freq <- 6
  ddt <- 0.001
  Ttot <- 2
  rate_mult <- 60
  nsubj <- 30

  simulate_amp <- function(amp) {
    oscore <- numeric(nsubj)
    fosc <- numeric(nsubj)
    exceeds_surrogate <- logical(nsubj)
    for (s in seq_len(nsubj)) {
      timeline <- seq(0, Ttot, by = ddt)
      ratetrace <- (1 + amp * sin(2 * pi * freq * timeline)) * rate_mult / length(timeline)
      events <- stats::runif(length(timeline)) < ratetrace
      RTs <- timeline[events]
      mct <- make_continuous_trace(RTs, dt = dt, warn = FALSE)
      os <- oscillation_score(mct$signal, fs = fs, flim = c(1, 20), warnings = FALSE)
      sur <- oscillation_score_surrogates(
        mct$signal, fs = fs, flim = c(1, 20), nrep = 20, fpeak = os$fosc, warnings = FALSE
      )
      oscore[s] <- os$oscore
      fosc[s] <- os$fosc
      exceeds_surrogate[s] <- is.finite(os$oscore) &&
        os$oscore > mean(sur$oscore_rp, na.rm = TRUE)
    }
    list(oscore = oscore, fosc = fosc, exceeds_surrogate = exceeds_surrogate)
  }

  flat <- simulate_amp(0)
  modulated <- simulate_amp(0.8)

  # 1. Oscillation score is substantially higher under rate modulation.
  expect_gt(
    mean(modulated$oscore, na.rm = TRUE),
    1.25 * mean(flat$oscore, na.rm = TRUE)
  )

  # 2. The recovered peak frequency clusters at the driving frequency.
  expect_lt(abs(stats::median(modulated$fosc, na.rm = TRUE) - freq), 1)
  expect_gte(sum(abs(modulated$fosc - freq) < 1, na.rm = TRUE), 20L)

  # 3. Most modulated subjects exceed their own surrogate distribution.
  expect_gte(sum(modulated$exceeds_surrogate), 20L)
})
