# Oscillation score with surrogate log-Z

Computes the oscillation score, generates a surrogate distribution, and
returns the log-Z score: \$\$Z = (log(O) - mean(log(O\_{rep}))) /
sd(log(O\_{rep}))\$\$.

## Usage

``` r
oscillation_score_z(
  signal,
  fs,
  flim,
  nrep = 500,
  alpha = 0.05,
  keep_trend = FALSE,
  trend_dist = c("gamma"),
  trend_ddt = NULL,
  trend_alpha = 0.05,
  surrogate_method = c("auto", "event_shuffle", "event_trend", "phase_randomized"),
  signal_mode = c("auto", "event", "continuous"),
  fcor = FALSE,
  taper = c("none", "hann", "hanning"),
  ci_nboot = 0,
  ci_level = 0.95,
  tidy = FALSE,
  ...
)

oscore_z(...)
```

## Arguments

- signal:

  numeric vector (binary or continuous).

- fs:

  sampling rate in Hz.

- flim:

  length-2 numeric vector `c(fmin, fmax)` in Hz.

- nrep:

  number of surrogates.

- alpha:

  significance level for one-tailed Z-thresholding.

- keep_trend, trend_dist, trend_ddt, trend_alpha:

  passed to surrogates.

- surrogate_method, signal_mode:

  passed to
  [`oscillation_score_surrogates`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_surrogates.md).

- fcor:

  logical; apply 1/f correction.

- taper:

  taper before FFT.

- ci_nboot:

  number of bootstrap draws for the log-Z confidence interval. Set to 0
  to skip CI estimation.

- ci_level:

  confidence level for bootstrap CI.

- tidy:

  logical; if TRUE return a one-row data.frame.

- ...:

  further arguments passed to
  [`oscillation_score`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score.md).

## Value

list (default) or data.frame if `tidy=TRUE`.
