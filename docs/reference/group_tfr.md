# Grouped time-frequency map via narrowband Hilbert decomposition

Computes a grouped time-frequency amplitude map by filtering each
subject series at a sequence of center frequencies and averaging the
resulting narrowband amplitudes.

## Usage

``` r
group_tfr(
  x,
  fs,
  freqs,
  bandwidth = NULL,
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

- freqs:

  numeric vector of center frequencies.

- bandwidth:

  optional half-bandwidth in Hz. If `NULL`, a default is derived from
  `freqs`.

- detrend:

  logical; if `TRUE`, detrend each subject series first.

- reflect:

  logical; if `TRUE`, reflect the series at both ends before filtering
  and trim back to the original span.

- edge:

  number of bins to trim from each edge after filtering.

- filtorder, demean, tol:

  passed to
  [`narrowband_hilbert`](https://bbuchsbaum.github.io/bosc/reference/narrowband_hilbert.md).

## Value

list with `observed` and `meta`.
