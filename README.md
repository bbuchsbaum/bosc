# bosc

R port of the [behavioral-oscillations](https://github.com/marijeterwal/behavioral-oscillations) MATLAB toolbox (Ter Wal et al.). Provides core functions for oscillation scoring, surrogate testing, phase analyses, PPC metrics, clustering, and plotting for behavioral and neural time series.

## Documentation

Full pkgdown site: <https://bbuchsbaum.github.io/bosc/>

### Vignettes

- [Get started with bosc](https://bbuchsbaum.github.io/bosc/articles/bosc.html) -- End-to-end introduction to the package workflow and core objects.
- [Oscillation Score Analysis](https://bbuchsbaum.github.io/bosc/articles/oscillation-score.html) -- Computing oscillation scores, surrogate testing, multi-band analysis, and narrowband Hilbert phase extraction.
- [PPC Analysis and Cluster Detection](https://bbuchsbaum.github.io/bosc/articles/ppc-clustering.html) -- Pairwise phase consistency, circular statistics, cluster-based permutation testing on time-frequency maps.
- [Workflow Wrappers](https://bbuchsbaum.github.io/bosc/articles/workflow-wrappers.html) -- High-level config-driven wrappers, batch analysis, phase-at-events pipeline, and tidy output helpers.

## Install

Install dependencies (`imager`, `signal`, `pracma`, `fitdistrplus`, `circular`, `ggplot2`, etc.) first, then:

```r
devtools::load_all("bosc")   # for development
# or install locally
devtools::install("bosc")
```

`imager` may require system image libraries (PNG/JPEG); install via your package manager before `install.packages("imager")`.

## Modules (current)

- Signal utils: `make_continuous_trace`, `autocorr_centered`, `spectral_peak`.
- O-score: `oscillation_score`, `oscillation_score_surrogates`.
- Workflow wrappers: `oscillation_score_config`, `oscillation_score_surrogates_config`, `oscillation_score_z`, `phase_at_events`.
- Tidy helpers: `oscore_tidy`, `oscore_spectrum`, `clusters_tidy`, `clusters_pixels`.
- Hilbert: `narrowband_hilbert`.
- PPC/Stats: `pairwise_phase_consistency`, `u_score_matrix`, `paired_tscore`, `nonparam_pval`.
- Clustering: `extract_clusters`, `detect_clusters` (imager-based watershed/extract, cluster splitting).
- Plotting: `plot_spectrum`, `plot_phase_hist`.

## Enhancements Relative to MATLAB Reference

This package closely follows the MATLAB reference implementation, with a few
small R-oriented enhancements for robustness, configurability, and workflow
ergonomics:

- `spectral_peak(fcor = TRUE)` uses an explicit, stable log-log fit/predict path
  in R for improved numerical stability.
- `oscillation_score()` supports explicit `signal_mode`:
  `event` uses event-rate-based effective bands; `continuous` uses
  variance-aware checks and Nyquist-limited upper bounds.
- `autocorr_centered()` now supports both direct and FFT-based linear
  autocorrelation; `method = "auto"` selects a fast FFT path for longer vectors
  while preserving centered `xcorr` semantics.
- `oscillation_score_surrogates()` supports explicit `surrogate_method`,
  including phase-randomized surrogates for continuous signals
  (amplitude-spectrum-preserving nulls), in addition to event shuffle/trend.
- `oscillation_score_z()` can return bootstrap confidence intervals for the
  log-Z statistic (`z_ci`) and exposes CI controls (`ci_nboot`, `ci_level`).
- `extract_clusters()` uses a direct 8-connectivity flood-fill implementation
  for consistent connected-component assignments.
- `circ_rayleigh()` uses weighted sample size (`sum(w)`) in weighted mode.
- `paired_tscore()` uses per-slice effective sample sizes under missing values.
- `detect_clusters(qval=...)` computes adjusted p-values, and optional
  `filter_by_q=TRUE` can filter clusters by adjusted significance.

## Testing

```r
devtools::test("bosc")
```

## References

This package is an R port of the MATLAB toolbox developed for:

**Ter Wal, M.**, Linde-Domingo, J., Lifanov, J., Roux, F., Kolibius, L., Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V., Chelvarajah, R., Staresina, B., Hanslmayr, S., & Wimber, M. (2021). Theta rhythmicity governs the timing of behavioural and hippocampal responses in humans specifically during memory-dependent tasks. *Nature Communications*, 12, 7048. <https://doi.org/10.1038/s41467-021-27323-3>

Dense-sampling trial-to-spectrum workflows in the grouped APIs and examples are also motivated by:

**Biba, T. M.**, Decker, A., Herrmann, B., Fukuda, K., Katz, C. N., Valiante, T. A., & Duncan, K. (2026). Episodic memory encoding fluctuates at a theta rhythm of 3-10 Hz. *Nature Human Behaviour*. <https://doi.org/10.1038/s41562-026-02416-5>

- **Original MATLAB code**: <https://github.com/marijeterwal/behavioral-oscillations>
- **Datasets**: <https://doi.org/10.6084/m9.figshare.c.5192567>
- **Video walkthrough** (CNS 2020 poster G183): <https://youtu.be/28kpbGDHLuo>

## Roadmap

- CI (R CMD check, lintr).
- Data I/O helpers for MAT/HDF5 and richer plotting wrappers.

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The palette family is provided via `params$family` (default 'lapis'). The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
