# bosc 0.1.0

* Initial CRAN release implementing behavioral oscillation scoring and phase
  analysis methods from Ter Wal et al. (2021), with grouped spectral and
  time-frequency workflows for trial-level behavioral data.
* Oscillation scores, surrogate testing, grouped spectra, phase TFRs, PPC metrics,
  clustering, and plotting utilities for behavioral and neural time series.

## Bug fixes and improvements

* `plot_simplebox()` now renders each group at its own position. Constant
  aesthetics had been placed inside `aes()`, so ggplot2's lazy evaluation drew
  every group at the last group's coordinates.
* `group_spectrum_test()`, `group_tfr_test()`, and the trial-level permutation
  workflow now report add-one (Phipson-Smyth) permutation p-values,
  `(1 + #{null >= obs}) / (1 + nrep)`, which are never exactly zero.
* `circ_mean()` and `circ_r()` drop non-finite angles for vector input, matching
  `circ_rayleigh()` so grouped phase summaries handle invalid subjects
  consistently.
* `aggregate_by_bin()` assembles its result without per-row `rbind`, materially
  speeding up the trial-level permutation workflows.
* `detect_clusters()` warns when a length-2 `smooth` requests anisotropic
  smoothing (only isotropic smoothing is supported). Documented the
  background-border convention used by `extract_clusters()` /
  `detect_clusters()` and the non-finite handling in the circular statistics.
