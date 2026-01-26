# Tidy spectrum from oscillation score or spectral peak output

Converts a result containing a spectrum into a long data.frame with
columns `freq` and `power`. Supports outputs from
[`oscillation_score`](oscillation_score.md) (fields `freqs`/`flimfft`)
and [`spectral_peak`](spectral_peak.md) (fields `fxx`/`spectrum`).

## Usage

``` r
oscore_spectrum(x)
```

## Arguments

- x:

  a spectrum-containing result list (or list of such results).

## Value

data.frame with columns `freq`, `power`, and optionally `id`.

## Details

If a list of results is provided, an `id` column is added.
