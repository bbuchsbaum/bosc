# Grouped spectrum with trial-level permutation nulls

Starts from a long trial table, aggregates within bins to form
per-subject series, and builds a null distribution by permuting bin
labels within each subject or subject-by-condition group before
re-aggregating.

## Usage

``` r
group_spectrum_test_trials(
  data,
  value,
  bin,
  subject,
  by = NULL,
  bins = NULL,
  fun = mean,
  na_rm = TRUE,
  complete = TRUE,
  incomplete = c("error", "drop"),
  fs,
  flim = NULL,
  null = c("shuffle_bins", "circular_shift_bins"),
  order_by = NULL,
  nrep = 1000,
  repair_incomplete = c("none", "resample"),
  max_resample = 100L,
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

- data:

  data.frame containing trial-wise values.

- value:

  character scalar naming the trial-level value column.

- bin:

  character scalar naming the bin column.

- subject:

  character scalar naming the subject column.

- by:

  optional character vector of additional grouping columns.

- bins:

  optional vector of target bins.

- fun:

  summary function used to aggregate trials within bins.

- na_rm:

  logical; if `TRUE`, remove missing trial values before aggregation.

- complete:

  logical; if `TRUE`, request explicit rows for missing group/bin
  combinations.

- incomplete:

  handling for incomplete observed grouped series; one of `"error"` or
  `"drop"`.

- fs:

  sampling rate in Hz along the bin axis.

- flim:

  optional length-2 frequency bounds.

- null:

  trial-level null generator; one of `"shuffle_bins"` or
  `"circular_shift_bins"`.

- order_by:

  optional within-group ordering column required for
  `"circular_shift_bins"`.

- nrep:

  number of null draws.

- repair_incomplete:

  how incomplete permuted groups are handled. `"none"` applies the same
  incomplete policy as the observed series. `"resample"` redraws each
  group's permutation until it yields a complete aggregated series.

- max_resample:

  maximum redraw attempts per group when
  `repair_incomplete = "resample"`.

- detrend:

  logical; if `TRUE`, detrend each series before spectral estimation.

- detrend_order:

  polynomial detrend order used when `detrend = TRUE`.

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

  optional random seed.

## Value

list with `observed`, `null`, `stats`, and `meta`. The `meta` field
includes the aggregated observed grouped-series object.

## Details

This workflow is intended for dense-sampling designs where trial timing
is varied across a regular grid and oscillatory structure is tested on
the resulting subject-by-bin series. The package implementation is
generic, but the memory-encoding application in Biba et al. (2026) is a
concrete example of this analysis pattern.

## References

Biba, T. M., Decker, A., Herrmann, B., Fukuda, K., Katz, C. N.,
Valiante, T. A., & Duncan, K. (2026). Episodic memory encoding
fluctuates at a theta rhythm of 3-10 Hz. *Nature Human Behaviour*.
[doi:10.1038/s41562-026-02416-5](https://doi.org/10.1038/s41562-026-02416-5)
