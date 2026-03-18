# Non-parametric p-values vs reference

Non-parametric p-values vs reference

## Usage

``` r
nonparam_pval(data, ref, alpha = 0.05)
```

## Arguments

- data:

  numeric vector.

- ref:

  numeric vector of reference samples.

- alpha:

  significance level.

## Value

list with `h` significance (0/1) and `p` p-values.

## Details

For each value in `data`, the p-value is computed as an empirical
one-sided tail proportion with respect to `ref`. Values below the
reference median use the lower tail (strictly smaller reference values),
while values at or above the median use the upper tail (strictly larger
reference values). Ties at the observed value are excluded from the tail
count by construction.
