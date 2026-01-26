# Paired t-score along a dimension

Paired t-score along a dimension

## Usage

``` r
paired_tscore(dat1, dat2 = NULL, dim = 1)
```

## Arguments

- dat1:

  array.

- dat2:

  optional array of same shape; if missing, one-sample t on `dat1`.

- dim:

  dimension to operate over.

## Value

list with `t` scores and `nu` degrees of freedom.
