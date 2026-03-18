# Build a grouped-series object from trial-level data

Aggregates trial values into regular bins for each subject or
subject-by-condition group, then reshapes the result into the
grouped-series structure accepted by
[`group_spectrum`](https://bbuchsbaum.github.io/bosc/reference/group_spectrum.md)
and related functions.

## Usage

``` r
make_grouped_series(
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
  out = NULL,
  meta = list()
)
```

## Arguments

- data:

  data.frame containing the trial-wise values.

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

  summary function applied within each bin.

- na_rm:

  logical; if `TRUE`, remove missing trial values before aggregation.

- complete:

  logical; if `TRUE`, request explicit rows for missing group/bin
  combinations.

- incomplete:

  handling for subject groups with missing or non-finite bin summaries.
  `"error"` stops, while `"drop"` removes incomplete groups.

- out:

  optional output measure name. Defaults to `value`.

- meta:

  optional named list of metadata to attach.

## Value

grouped-series list with `data`, `subjects`, `bins`, `measure`, `meta`,
plus `aggregated` and `group_keys`.
