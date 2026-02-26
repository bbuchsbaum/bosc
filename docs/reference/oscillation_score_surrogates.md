# Surrogate oscillation scores (shuffle or trend-preserving)

Generates surrogate signals and computes oscillation scores for
reference distributions. Fixes edge cases from the MATLAB implementation
(e.g., negative p-values, missing sums).

## Usage

``` r
oscillation_score_surrogates(
  signal,
  fs,
  flim,
  nrep,
  fpeak,
  keep_trend = FALSE,
  trend_dist = c("gamma"),
  trend_ddt = NULL,
  trend_alpha = 0.05,
  surrogate_method = c("auto", "event_shuffle", "event_trend", "phase_randomized"),
  signal_mode = c("auto", "event", "continuous"),
  warnings = TRUE,
  fcor = FALSE,
  taper = c("none", "hann", "hanning")
)
```

## Arguments

- signal:

  numeric vector.

- fs:

  sampling rate in Hz.

- flim:

  length-2 numeric vector `c(fmin, fmax)`.

- nrep:

  number of surrogates.

- fpeak:

  peak frequency of original signal (Hz); required for windowing.

- keep_trend:

  logical; if TRUE fit trend distribution instead of shuffling.

- trend_dist:

  character vector of distributions for
  [`fitdistrplus::fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html).

- trend_ddt:

  step (seconds) for trend generation; defaults to `0.5/fs`.

- trend_alpha:

  p-value threshold to accept fit.

- surrogate_method:

  surrogate generator; one of `"auto"`, `"event_shuffle"`,
  `"event_trend"`, or `"phase_randomized"`. In `"auto"` mode, continuous
  signals default to phase randomization, while event-like signals use
  shuffle or trend fitting depending on `keep_trend`.

- signal_mode:

  one of `"auto"`, `"event"`, or `"continuous"`.

- warnings:

  logical; emit warnings.

- fcor:

  logical; apply 1/f correction in surrogate peak estimation.

- taper:

  taper applied before FFT; one of `"none"`, `"hann"`, `"hanning"`.

## Value

list with `oscore_rp`, `fosc_rp`, `trendfit`, `surrogate_method`,
`signal_mode`, and (optionally) `signrep`.
