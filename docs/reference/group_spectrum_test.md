# Grouped spectrum with resampled null testing

Builds a null distribution for the group-average spectrum by resampling
each subject series, then computes empirical p-values, z-scores, and
optional FDR correction across frequencies.

## Usage

``` r
group_spectrum_test(
  x,
  fs,
  flim = NULL,
  null = c("shuffle_labels", "circular_shift"),
  nrep = 1000,
  detrend = TRUE,
  detrend_order = 1L,
  pad_to = NULL,
  spectrum = c("amplitude", "raw_power"),
  taper = c("none", "hann", "hanning"),
  p_adjust = c("none", "fdr"),
  fcor = FALSE,
  seed = NULL
)
```

## Arguments

- x:

  numeric matrix with shape `subjects x bins`, or a grouped-series
  object.

- fs:

  sampling rate in Hz along the bin axis.

- flim:

  optional length-2 frequency bounds.

- null:

  null generator; one of `"shuffle_labels"` or `"circular_shift"`.

- nrep:

  number of null draws.

- detrend:

  logical; if `TRUE`, remove a linear trend from each subject series
  before spectral estimation.

- detrend_order:

  non-negative integer polynomial order used when `detrend = TRUE`.

- pad_to:

  optional FFT length used for zero-padding each subject series before
  spectral estimation.

- spectrum:

  spectrum scaling passed to
  [`spectral_peak`](https://bbuchsbaum.github.io/bosc/reference/spectral_peak.md).

- taper:

  taper passed to
  [`spectral_peak`](https://bbuchsbaum.github.io/bosc/reference/spectral_peak.md).

- p_adjust:

  p-value adjustment; one of `"none"` or `"fdr"`.

- fcor:

  logical; apply 1/f correction to each subject spectrum.

- seed:

  optional random seed for reproducible null draws.

## Value

list with `observed`, `null`, `stats`, and `meta`.
