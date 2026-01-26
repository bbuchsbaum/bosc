# Connected component labeling for thresholded data

Thresholds 2D data and labels connected components (8-connectivity).

## Usage

``` r
extract_clusters(data, thres = NULL)
```

## Arguments

- data:

  matrix of values.

- thres:

  numeric threshold; if NULL, assumes binary input.

## Value

integer matrix of labels (0 for below threshold).
