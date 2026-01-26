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
  fcor = FALSE,
  taper = c("none", "hann", "hanning"),
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

- fcor:

  logical; apply 1/f correction.

- taper:

  taper before FFT.

- tidy:

  logical; if TRUE return a one-row data.frame.

- ...:

  further arguments passed to
  [`oscillation_score`](oscillation_score.md).

## Value

list (default) or data.frame if `tidy=TRUE`.
