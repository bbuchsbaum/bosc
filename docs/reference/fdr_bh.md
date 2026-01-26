# FDR Correction (Benjamini-Hochberg)

Wrapper around
[`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for FDR
correction using the Benjamini-Hochberg procedure. Provides interface
similar to MATLAB `fdr_bh`.

## Usage

``` r
fdr_bh(pvals, q = 0.05, method = "BH")
```

## Arguments

- pvals:

  numeric vector of p-values.

- q:

  desired false discovery rate (default 0.05).

- method:

  adjustment method passed to `p.adjust`. Default "BH".

## Value

list with:

- h:

  logical vector indicating significant tests after FDR correction

- crit_p:

  critical p-value threshold (largest p-value deemed significant)

- adj_p:

  adjusted p-values

## References

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. Journal of the
Royal Statistical Society B, 57, 289-300.

## Examples

``` r
pvals <- c(0.001, 0.01, 0.03, 0.04, 0.05, 0.1, 0.5)
result <- fdr_bh(pvals, q = 0.05)
result$h
#> [1]  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE
```
