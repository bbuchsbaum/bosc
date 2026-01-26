# Mean Resultant Vector Length

Computes the mean resultant vector length for circular data, a measure
of concentration (0 = uniform, 1 = perfectly aligned). Port of MATLAB
`circ_r` from the Circular Statistics Toolbox (Berens, 2009).

## Usage

``` r
circ_r(alpha, w = NULL, d = 0, dim = 1)
```

## Arguments

- alpha:

  numeric vector or matrix of angles in radians.

- w:

  optional weights (same dimensions as `alpha`). If NULL, uniform
  weights.

- d:

  optional spacing of bin centers for binned data (radians). If
  provided, applies correction factor for bias in estimation of r.

- dim:

  integer; dimension along which to compute (default 1).

## Value

Mean resultant length (scalar for vector input, vector/array otherwise).

## References

Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics.
Journal of Statistical Software, 31(10), 1-21.

Zar JH (2010). Biostatistical Analysis. 5th ed. Prentice Hall. (Equation
26.16)

## See also

[`circ_mean`](circ_mean.md), [`circ_vtest`](circ_vtest.md)

## Examples

``` r
# Highly concentrated data
angles <- c(0, 0.1, -0.1, 0.05)
circ_r(angles)  # Should be close to 1
#> [1] 0.9972679

# Uniformly distributed data
angles_uniform <- seq(0, 2*pi, length.out = 100)
circ_r(angles_uniform)  # Should be close to 0
#> [1] 0.01
```
