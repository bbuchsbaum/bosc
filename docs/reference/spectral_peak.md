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
  pad_to = NULL,
  spectrum = c("amplitude", "raw_power"),
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

- pad_to:

  optional integer length for zero-padding before the FFT. When `NULL`,
  no zero-padding is applied.

- spectrum:

  spectrum scaling. `"amplitude"` returns the current normalized
  single-sided amplitude spectrum. `"raw_power"` returns raw FFT power
  `|FFT(x)|^2` on the retained positive-frequency bins.

- taper:

  taper to apply; one of `"none"`, `"hann"`, `"hanning"`.

## Value

list with `freq` (peak frequency or NA), `fxx` (freq axis), and
`spectrum` (power values aligned with `fxx`).
