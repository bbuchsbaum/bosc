# Mann-Whitney U Z-score per element vs. reference repeats

Mann-Whitney U Z-score per element vs. reference repeats

## Usage

``` r
u_score_matrix(dat, ref, alpha = 0.05)
```

## Arguments

- dat:

  array of data.

- ref:

  array where first dimension is repetitions.

- alpha:

  significance threshold.

## Value

list with `Z` (same shape as `dat`) and `sgnf` (significance mask).
