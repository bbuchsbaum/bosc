# Grouped time-frequency map with resampled null testing

Computes an empirical null distribution for grouped narrowband amplitude
maps and returns cell-wise statistics with optional cluster summaries.

## Usage

``` r
group_tfr_test(
  x,
  fs,
  freqs,
  bandwidth = NULL,
  null = c("shuffle_labels", "circular_shift"),
  nrep = 1000,
  detrend = TRUE,
  reflect = FALSE,
  edge = 0,
  filtorder = 2,
  demean = TRUE,
  tol = 100,
  p_adjust = c("none", "fdr"),
  cluster = TRUE,
  cluster_threshold = stats::qnorm(0.975),
  seed = NULL
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

  optional half-bandwidth in Hz.

- null:

  null generator; one of `"shuffle_labels"` or `"circular_shift"`.

- nrep:

  number of null draws.

- detrend, reflect, edge, filtorder, demean, tol:

  passed to
  [`group_tfr`](https://bbuchsbaum.github.io/bosc/reference/group_tfr.md).

- p_adjust:

  p-value adjustment; one of `"none"` or `"fdr"`.

- cluster:

  logical; if `TRUE`, detect suprathreshold clusters in the z-map.

- cluster_threshold:

  z-threshold used for cluster detection.

- seed:

  optional random seed.

## Value

list with `observed`, `null`, `stats`, and `meta`.
