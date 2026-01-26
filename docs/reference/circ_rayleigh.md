# Rayleigh Test for Non-uniformity

Tests whether a sample of circular data is uniformly distributed. This
is a wrapper around the resultant length statistic.

## Usage

``` r
circ_rayleigh(alpha, w = NULL)
```

## Arguments

- alpha:

  numeric vector of angles in radians.

- w:

  optional weights (same length as `alpha`).

## Value

list with:

- pval:

  p-value of the Rayleigh test

- z:

  Rayleigh's Z statistic

## References

Zar JH (2010). Biostatistical Analysis. 5th ed. Prentice Hall.

## Examples

``` r
# Uniform data - expect high p-value
uniform_angles <- runif(100, 0, 2*pi)
circ_rayleigh(uniform_angles)
#> $pval
#> [1] 0.1835458
#> 
#> $z
#> [1] 1.696592
#> 

# Concentrated data - expect low p-value
concentrated <- rnorm(100, mean = 0, sd = 0.3)
circ_rayleigh(concentrated)
#> $pval
#> [1] 3.258676e-38
#> 
#> $z
#> [1] 91.50884
#> 
```
