# Configuration-list wrapper for surrogate O-scores

Companion to [`oscillation_score_config`](oscillation_score_config.md)
for surrogate testing. Accepts a config list (named like
[`oscillation_score_surrogates`](oscillation_score_surrogates.md)) and
forwards it to the surrogate function.

## Usage

``` r
oscillation_score_surrogates_config(config, signal)

oscillation_score_stats(cfg, signal)
```

## Arguments

- config:

  list of settings.

- signal:

  numeric vector.

## Value

list as returned by
[`oscillation_score_surrogates`](oscillation_score_surrogates.md).

## Details

Recognized fields include `fs`, `flim`, `nrep`, `fpeak`, `keep_trend`,
`trend_dist`, `trend_ddt`, `trend_alpha`, `warnings`, `fcor`, and
`taper`.
