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

- method:

  autocorrelation backend: `"direct"` (explicit summation), `"fft"`
  (FFT-based linear autocorrelation), or `"auto"` (uses `"direct"` for
  short vectors and `"fft"` for longer vectors).

## Value

numeric vector of length `2*length(x) - 1`.
