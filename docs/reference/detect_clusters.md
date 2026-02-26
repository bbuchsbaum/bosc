# Detect clusters with optional splitting

Detect clusters with optional splitting

## Usage

``` r
detect_clusters(
  data,
  threshold,
  method = c("watershed", "extract"),
  smooth = NULL,
  split_clusters = FALSE,
  peak_height = seq(0.05, 0.5, by = 0.05),
  min_new_cluster_size = 0.1,
  freqweight = 1,
  timeweight = 1,
  freq_axis = NULL,
  time_axis = NULL,
  data2 = NULL,
  dataRef = NULL,
  qval = NULL,
  filter_by_q = FALSE
)
```

## Arguments

- data:

  matrix of primary values.

- threshold:

  numeric; binarization threshold.

- method:

  "extract" or "watershed".

- smooth:

  optional numeric or length-2 vector for Gaussian smoothing (pixels).

- split_clusters:

  logical; attempt to split multi-peak clusters.

- peak_height:

  numeric vector of relative peak heights for splitting.

- min_new_cluster_size:

  minimum fraction of original cluster to keep split.

- freqweight, timeweight:

  weights for distance reassignment.

- freq_axis, time_axis:

  vectors for frequencies/times corresponding to rows/cols.

- data2:

  optional secondary matrix for peak value extraction.

- dataRef:

  optional reference matrix for computing reference statistics per
  cluster.

- qval:

  optional FDR q for p-value adjustments.

- filter_by_q:

  logical; if TRUE and `qval` is provided, keep only clusters with at
  least one pixel passing `pAdj <= qval`.

## Value

list of cluster structs with coordinates and stats.
