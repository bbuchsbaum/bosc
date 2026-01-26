# Tidy per-pixel cluster data

Returns a long data.frame with one row per cluster pixel. Useful for
plotting or post-hoc summaries.

## Usage

``` r
clusters_pixels(clusters)
```

## Arguments

- clusters:

  list of clusters as returned by
  [`detect_clusters`](detect_clusters.md).

## Value

data.frame with columns `cluster`, `time`, `freq`, `Zscore`, and `pAdj`
(if present).
