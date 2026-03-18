# Permute trial labels within groups

Shuffles or circularly shifts a label column within groups. This is
useful for trial-level permutation nulls where the trial values are left
unchanged but their bin assignments are resampled.

## Usage

``` r
permute_trial_labels(
  data,
  column,
  by = NULL,
  method = c("shuffle", "circular_shift"),
  order_by = NULL,
  out = NULL,
  seed = NULL
)
```

## Arguments

- data:

  data.frame containing the label column.

- column:

  character scalar naming the column to permute.

- by:

  optional character vector of grouping columns.

- method:

  permutation method. `"shuffle"` samples without replacement.
  `"circular_shift"` preserves local order after sorting by `order_by`.

- order_by:

  optional character scalar naming the within-group ordering column
  required for `"circular_shift"`.

- out:

  optional output column name. Defaults to `column`.

- seed:

  optional random seed.

## Value

data.frame with the permuted label column.
