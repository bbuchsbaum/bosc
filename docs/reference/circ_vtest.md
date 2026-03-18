# V-test for Non-uniformity with Specified Mean Direction

Computes the V-test for circular data, testing whether the population is
non-uniformly distributed with a specified mean direction. Port of
MATLAB `circ_vtest` from the Circular Statistics Toolbox (Berens, 2009).

## Usage

``` r
circ_vtest(alpha, dir, w = NULL, d = 0)
```

## Arguments

- alpha:

  numeric vector of angles in radians.

- dir:

  suspected mean direction in radians.

- w:

  optional weights for binned data (same length as `alpha`).

- d:

  optional bin spacing for binned data (radians) for bias correction.

## Value

list with:

- pval:

  p-value of the V-test (one-tailed)

- v:

  V statistic

## Details

The V-test is more powerful than the Rayleigh test when there is reason
to believe in a specific mean direction.

## References

Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics.
Journal of Statistical Software, 31(10), 1-21.

Zar JH (2010). Biostatistical Analysis. 5th ed. Prentice Hall.

## See also

[`circ_mean`](https://bbuchsbaum.github.io/bosc/reference/circ_mean.md),
[`circ_r`](https://bbuchsbaum.github.io/bosc/reference/circ_r.md)

## Examples

``` r
# Test if data clusters around 0
angles <- rnorm(50, mean = 0, sd = 0.5)
result <- circ_vtest(angles, dir = 0)
result$pval
#> [1] 0
```
