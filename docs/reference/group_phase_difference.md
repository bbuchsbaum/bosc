# Grouped phase difference over bins

Computes per-bin phase differences between two matched grouped series
objects or matrices, then summarizes phase concentration of those
differences across subjects.

## Usage

``` r
group_phase_difference(
  x1,
  x2,
  fs,
  freq,
  bandwidth = 0.5,
  mode = c("signed_pi", "wrapped_2pi"),
  detrend = TRUE,
  reflect = FALSE,
  edge = 0,
  filtorder = 2,
  demean = TRUE,
  tol = 100
)
```

## Arguments

- x1, x2:

  numeric matrices with identical shape, or grouped-series objects with
  matching `subjects` and `bins`.

- fs:

  sampling rate in Hz along the bin axis.

- freq:

  center frequency in Hz.

- bandwidth:

  half-bandwidth in Hz.

- mode:

  one of `"signed_pi"` or `"wrapped_2pi"`.

- detrend, reflect, edge, filtorder, demean, tol:

  preprocessing controls.

## Value

list with `observed`, `stats`, and `meta`.
