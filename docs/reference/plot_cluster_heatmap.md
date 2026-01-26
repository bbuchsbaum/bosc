# Plot Cluster Heatmap

Visualize detected clusters as a heatmap with optional cluster boundary
overlays.

## Usage

``` r
plot_cluster_heatmap(
  data,
  clusters = NULL,
  freq_axis = NULL,
  time_axis = NULL,
  zlim = NULL,
  palette = "RdBu"
)
```

## Arguments

- data:

  matrix of values (e.g., Z-scores).

- clusters:

  optional list of cluster objects from `detect_clusters`.

- freq_axis:

  optional frequency axis for y-axis labels.

- time_axis:

  optional time axis for x-axis labels.

- zlim:

  optional length-2 vector for color scale limits.

- palette:

  color palette name (default "RdBu").

## Value

ggplot object.

## Examples

``` r
mat <- matrix(rnorm(400), nrow = 20)
mat[8:12, 8:12] <- mat[8:12, 8:12] + 3
plot_cluster_heatmap(mat)
```
