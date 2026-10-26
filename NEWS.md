# bosc 0.1.0

* Initial CRAN release implementing behavioral oscillation scoring and phase
  analysis methods from Ter Wal et al. (2021), with grouped spectral and
  time-frequency workflows for trial-level behavioral data.
* Oscillation scores, surrogate testing, grouped spectra, phase TFRs, PPC metrics,
  clustering, and plotting utilities for behavioral and neural time series.
* `make_continuous_trace()` now removes floating-point roundoff (including tiny
  negative values) from the convolved trace, so event trains stay non-negative
  and sparse. Previously `signal_mode = "auto"` could misclassify such traces
  as continuous, which made `oscillation_score_surrogates()` use
  phase-randomised surrogates that preserve the autocorrelation and hence the
  oscillation score. Event-like detection now also tolerates roundoff.
* Added examples to the main exported functions.
