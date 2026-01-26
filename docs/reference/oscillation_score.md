# Oscillation Score for binary or continuous data

Port of the MATLAB `oscillationScore` algorithm with robustness fixes.
Computes an autocorrelogram, removes the central peak, finds the
spectral peak, and returns the oscillation score.

## Usage

``` r
oscillation_score(
  signal,
  fs,
  flim,
  quantlim = NULL,
  smoothach = TRUE,
  smoothwind = NULL,
  peakwind = NULL,
  thresangle = 10,
  mincycles = 3,
  minfreqbandwidth = NULL,
  fpeak = NULL,
  warnings = TRUE,
  fcor = FALSE,
  taper = c("none", "hann", "hanning"),
  plot = FALSE
)
```

## Arguments

- signal:

  numeric vector (binary or continuous).

- fs:

  sampling rate in Hz.

- flim:

  length-2 numeric vector `c(fmin, fmax)` in Hz.

- quantlim:

  optional length-2 quantile bounds for trimming nonzero indices.

- smoothach:

  logical; smooth the ACH with a Gaussian kernel.

- smoothwind:

  smoothing kernel width (seconds); default matches MATLAB heuristic.

- peakwind:

  slow-kernel width for peak removal (seconds); default matches MATLAB
  heuristic.

- thresangle:

  angular threshold for peak detection (degrees).

- mincycles:

  minimum cycles to include at `fmin`.

- minfreqbandwidth:

  optional minimum bandwidth; returns NA if unmet.

- fpeak:

  optional frequency of interest; overrides peak search.

- warnings:

  logical; emit warnings.

- fcor:

  logical; apply 1/f correction in spectral peak estimation.

- taper:

  taper applied before FFT; one of `"none"`, `"hann"`, `"hanning"`.

- plot:

  logical; if TRUE, plots intermediate spectra (optional, requires
  ggplot2).

## Value

list with `oscore`, `fosc`, `flim`, `flimfft`, `freqs`.
