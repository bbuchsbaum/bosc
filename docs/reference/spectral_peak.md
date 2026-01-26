# Find spectral peak within bounds

Computes the FFT power spectrum (single-sided), optional 1/f correction,
and returns the frequency of the highest peak within bounds.

## Usage

``` r
spectral_peak(
  data,
  fs,
  flim = NULL,
  fcor = FALSE,
  taper = c("none", "hann", "hanning")
)
```

## Arguments

- data:

  numeric vector.

- fs:

  sampling rate (Hz).

- flim:

  optional length-2 frequency bounds `c(fmin, fmax)` in Hz.

- fcor:

  logical; apply crude 1/f correction.

- taper:

  taper to apply; one of `"none"`, `"hann"`, `"hanning"`.

## Value

list with `freq` (peak frequency or NA), `fxx` (freq axis), and
`spectrum` (power values aligned with `fxx`).
