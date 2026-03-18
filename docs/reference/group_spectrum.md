# Grouped spectrum for binned subject-level series

Computes per-subject spectra and their group average for a regular
matrix of subject-by-bin values or a grouped-series object.

## Usage

``` r
group_spectrum(
  x,
  fs,
  flim = NULL,
  detrend = TRUE,
  detrend_order = 1L,
  pad_to = NULL,
  spectrum = c("amplitude", "raw_power"),
  taper = c("none", "hann", "hanning"),
  average = TRUE,
  fcor = FALSE
)
```

## Arguments

- x:

  numeric matrix with shape `subjects x bins`, or a list with fields
  `data`, `subjects`, and `bins`.

- fs:

  sampling rate in Hz along the bin axis.

- flim:

  optional length-2 frequency bounds.

- detrend:

  logical; if `TRUE`, remove a linear trend from each subject series
  before spectral estimation.

- detrend_order:

  non-negative integer polynomial order used when `detrend = TRUE`. The
  default `1` removes a linear trend, while higher orders remove
  higher-order polynomial structure.

- pad_to:

  optional FFT length used for zero-padding each subject series before
  spectral estimation.

- spectrum:

  spectrum scaling passed to
  [`spectral_peak`](https://bbuchsbaum.github.io/bosc/reference/spectral_peak.md).

- taper:

  taper passed to
  [`spectral_peak`](https://bbuchsbaum.github.io/bosc/reference/spectral_peak.md).

- average:

  logical; if `TRUE`, include the group-average spectrum.

- fcor:

  logical; apply 1/f correction to each subject spectrum.

## Value

list with `observed` and `meta` fields.
