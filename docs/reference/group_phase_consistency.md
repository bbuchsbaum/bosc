# Grouped phase consistency over bins

Computes per-bin phase concentration across subjects after narrowband
Hilbert decomposition of subject-level series.

## Usage

``` r
group_phase_consistency(
  x,
  fs,
  freq,
  bandwidth = 0.5,
  detrend = TRUE,
  reflect = FALSE,
  edge = 0,
  filtorder = 2,
  demean = TRUE,
  tol = 100
)
```

## Arguments

- x:

  numeric matrix with shape `subjects x bins`, or a grouped-series
  object.

- fs:

  sampling rate in Hz along the bin axis.

- freq:

  center frequency in Hz.

- bandwidth:

  half-bandwidth in Hz.

- detrend, reflect, edge, filtorder, demean, tol:

  preprocessing controls.

## Value

list with `observed`, `stats`, and `meta`.
