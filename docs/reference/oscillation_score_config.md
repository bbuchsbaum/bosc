# Configuration-list wrapper for oscillation scores

Accepts a list of parameters (named like
[`oscillation_score`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score.md))
and forwards them to the core O-score function. Useful when your
analysis builds configs programmatically.

## Usage

``` r
oscillation_score_config(config, signal)

oscillation_score_cfg(cfg, signal)
```

## Arguments

- config:

  list of settings.

- signal:

  numeric vector (binary or continuous).

- cfg:

  list of settings (deprecated alias for `config`).

## Value

list as returned by
[`oscillation_score`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score.md).

## Details

Recognized fields include `fs`, `flim`, `quantlim`, `smoothach`,
`smoothwind`, `peakwind`, `thresangle`, `mincycles`, `minfreqbandwidth`,
`fpeak`, `warnings`, `fcor`, `taper`, `signal_mode`, and `plot`.
