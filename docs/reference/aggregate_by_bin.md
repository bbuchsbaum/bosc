# Aggregate trial values into regularized bins

Aggregates trial-level values within bins and optional grouping
variables, producing a tidy table suitable for grouped time-series
analysis.

## Usage

``` r
aggregate_by_bin(
  data,
  value,
  bin,
  by = NULL,
  bins = NULL,
  fun = mean,
  na_rm = TRUE,
  complete = TRUE,
  out = NULL
)
```

## Arguments

- data:

  data.frame containing the value and bin columns.

- value:

  character scalar naming the value column to aggregate.

- bin:

  character scalar naming the bin column.

- by:

  optional character vector of grouping columns.

- bins:

  optional vector of target bins. If `NULL`, uses the sorted unique
  non-missing values observed in `data[[bin]]`.

- fun:

  summary function applied within each bin.

- na_rm:

  logical; if `TRUE`, remove missing values before applying `fun`.

- complete:

  logical; if `TRUE`, return explicit rows for missing group/bin
  combinations with `NA` summaries.

- out:

  optional output column name. Defaults to `value`.

## Value

data.frame with grouping columns, `bin`, the aggregated value, and count
columns `n` and `n_nonmissing`.
