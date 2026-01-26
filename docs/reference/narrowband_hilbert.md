# Narrowband filtering and Hilbert transform

Filters data in a narrow frequency band and computes the analytic signal
via the Hilbert transform. Mirrors MATLAB `narrowbandHilbert` with
stability checks and optional demeaning.

## Usage

``` r
narrowband_hilbert(data, fs, freqlim, tol = 100, filtorder = 2, demean = TRUE)
```

## Arguments

- data:

  numeric vector.

- fs:

  sampling rate (Hz).

- freqlim:

  length-2 numeric vector `c(fmin, fmax)` in Hz.

- tol:

  tolerance for filter stability (max absolute value before flagging).

- filtorder:

  filter order (integer).

- demean:

  logical; remove mean before filtering.

## Value

list with `analytic` (complex), `filtered` (numeric).
