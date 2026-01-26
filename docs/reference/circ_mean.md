# Circular Mean Direction

Computes the mean direction for circular data. Port of MATLAB
`circ_mean` from the Circular Statistics Toolbox (Berens, 2009).

## Usage

``` r
circ_mean(alpha, w = NULL, dim = 1)
```

## Arguments

- alpha:

  numeric vector or matrix of angles in radians.

- w:

  optional weights (same dimensions as `alpha`). If NULL, uniform
  weights.

- dim:

  integer; dimension along which to compute mean (default 1).

## Value

If `alpha` is a vector, returns a single mean direction. If matrix,
returns a vector of means computed along `dim`.

## References

Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics.
Journal of Statistical Software, 31(10), 1-21.

Fisher NI (1993). Statistical Analysis of Circular Data. Cambridge
University Press.

## See also

[`circ_r`](circ_r.md), [`circ_vtest`](circ_vtest.md)

## Examples

``` r
angles <- c(0, pi/4, pi/2)
circ_mean(angles)
#> [1] 0.7853982

# Weighted mean
circ_mean(angles, w = c(1, 2, 1))
#> [1] 0.7853982
```
