# Tidy summary of detected clusters

Converts the list returned by
[`detect_clusters`](https://bbuchsbaum.github.io/bosc/reference/detect_clusters.md)
into a data.frame with one row per cluster.

## Usage

``` r
clusters_tidy(clusters)
```

## Arguments

- clusters:

  list of clusters as returned by
  [`detect_clusters`](https://bbuchsbaum.github.io/bosc/reference/detect_clusters.md).

## Value

data.frame with cluster summaries.
