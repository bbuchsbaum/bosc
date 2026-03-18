# Tidy summary of oscillation score results

Converts an
[`oscillation_score`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score.md)
or
[`oscillation_score_z`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_z.md)
result into a one-row data.frame. If a list of results is provided,
returns one row per element.

## Usage

``` r
oscore_tidy(x)
```

## Arguments

- x:

  a result list (or list of result lists).

## Value

data.frame with summary columns.
