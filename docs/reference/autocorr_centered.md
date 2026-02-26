# Centered autocorrelation (unnormalized)

Computes an unnormalized autocorrelation with zero lag at the center,
equivalent to MATLAB `xcorr(x)`.

## Usage

``` r
autocorr_centered(x, method = c("auto", "direct", "fft"))
```

## Arguments

- x:

  numeric vector.

## Value

numeric vector of length `2*length(x) - 1`.
