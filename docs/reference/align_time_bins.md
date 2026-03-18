# Align observed times to a target bin grid

Maps observed timing values to the nearest target bin. This is useful
when nominal bins are known in advance but measured times contain small
timing deviations.

## Usage

``` r
align_time_bins(
  data,
  time,
  bins,
  tolerance = NULL,
  method = c("nearest"),
  out = NULL
)
```

## Arguments

- data:

  data.frame containing the timing column.

- time:

  character scalar naming the observed time column.

- bins:

  numeric vector of target bin values.

- tolerance:

  optional non-negative scalar. If provided, values farther than
  `tolerance` from the nearest target bin are assigned `NA`.

- method:

  matching method. Currently only `"nearest"` is supported.

- out:

  optional output column name. Defaults to `paste0(time, "_bin")`.

## Value

data.frame with an added aligned-bin column.
