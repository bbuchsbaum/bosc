# Package index

## Signal and O-score

- [`make_continuous_trace()`](https://bbuchsbaum.github.io/bosc/reference/make_continuous_trace.md)
  : Make a continuous trace from event times
- [`autocorr_centered()`](https://bbuchsbaum.github.io/bosc/reference/autocorr_centered.md)
  : Centered autocorrelation (unnormalized)
- [`spectral_peak()`](https://bbuchsbaum.github.io/bosc/reference/spectral_peak.md)
  : Find spectral peak within bounds
- [`oscillation_score()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score.md)
  : Oscillation Score for binary or continuous data
- [`oscillation_score_surrogates()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_surrogates.md)
  : Surrogate oscillation scores (shuffle or trend-preserving)

## Workflow Wrappers

- [`oscillation_score_z()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_z.md)
  [`oscore_z()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_z.md)
  : Oscillation score with surrogate log-Z
- [`oscillation_score_config()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_config.md)
  [`oscillation_score_cfg()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_config.md)
  : Configuration-list wrapper for oscillation scores
- [`oscillation_score_surrogates_config()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_surrogates_config.md)
  [`oscillation_score_stats()`](https://bbuchsbaum.github.io/bosc/reference/oscillation_score_surrogates_config.md)
  : Configuration-list wrapper for surrogate O-scores
- [`phase_at_events()`](https://bbuchsbaum.github.io/bosc/reference/phase_at_events.md)
  : Extract narrowband phase at event times

## Tidy Helpers

- [`oscore_tidy()`](https://bbuchsbaum.github.io/bosc/reference/oscore_tidy.md)
  : Tidy summary of oscillation score results
- [`oscore_spectrum()`](https://bbuchsbaum.github.io/bosc/reference/oscore_spectrum.md)
  : Tidy spectrum from oscillation score or spectral peak output
- [`clusters_tidy()`](https://bbuchsbaum.github.io/bosc/reference/clusters_tidy.md)
  : Tidy summary of detected clusters
- [`clusters_pixels()`](https://bbuchsbaum.github.io/bosc/reference/clusters_pixels.md)
  : Tidy per-pixel cluster data

## Phase and Hilbert

- [`narrowband_hilbert()`](https://bbuchsbaum.github.io/bosc/reference/narrowband_hilbert.md)
  : Narrowband filtering and Hilbert transform

## Circular Statistics

- [`circ_mean()`](https://bbuchsbaum.github.io/bosc/reference/circ_mean.md)
  : Circular Mean Direction
- [`circ_r()`](https://bbuchsbaum.github.io/bosc/reference/circ_r.md) :
  Mean Resultant Vector Length
- [`circ_rayleigh()`](https://bbuchsbaum.github.io/bosc/reference/circ_rayleigh.md)
  : Rayleigh Test for Non-uniformity
- [`circ_vtest()`](https://bbuchsbaum.github.io/bosc/reference/circ_vtest.md)
  : V-test for Non-uniformity with Specified Mean Direction

## PPC and Statistical Testing

- [`pairwise_phase_consistency()`](https://bbuchsbaum.github.io/bosc/reference/pairwise_phase_consistency.md)
  : Pairwise Phase Consistency (1D)
- [`u_score_matrix()`](https://bbuchsbaum.github.io/bosc/reference/u_score_matrix.md)
  : Mann-Whitney U Z-score per element vs. reference repeats
- [`paired_tscore()`](https://bbuchsbaum.github.io/bosc/reference/paired_tscore.md)
  : Paired t-score along a dimension
- [`nonparam_pval()`](https://bbuchsbaum.github.io/bosc/reference/nonparam_pval.md)
  : Non-parametric p-values vs reference
- [`fdr_bh()`](https://bbuchsbaum.github.io/bosc/reference/fdr_bh.md) :
  FDR Correction (Benjamini-Hochberg)

## Clustering

- [`extract_clusters()`](https://bbuchsbaum.github.io/bosc/reference/extract_clusters.md)
  : Connected component labeling for thresholded data
- [`detect_clusters()`](https://bbuchsbaum.github.io/bosc/reference/detect_clusters.md)
  : Detect clusters with optional splitting

## Grouped Data Preparation

- [`align_time_bins()`](https://bbuchsbaum.github.io/bosc/reference/align_time_bins.md)
  : Align observed times to a target bin grid
- [`residualize_trials()`](https://bbuchsbaum.github.io/bosc/reference/residualize_trials.md)
  : Residualize trial-wise responses against nuisance terms
- [`aggregate_by_bin()`](https://bbuchsbaum.github.io/bosc/reference/aggregate_by_bin.md)
  : Aggregate trial values into regularized bins

## Grouped Analysis

- [`group_spectrum()`](https://bbuchsbaum.github.io/bosc/reference/group_spectrum.md)
  : Grouped spectrum for binned subject-level series
- [`group_spectrum_test()`](https://bbuchsbaum.github.io/bosc/reference/group_spectrum_test.md)
  : Grouped spectrum with resampled null testing
- [`group_tfr()`](https://bbuchsbaum.github.io/bosc/reference/group_tfr.md)
  : Grouped time-frequency map via narrowband Hilbert decomposition
- [`group_tfr_test()`](https://bbuchsbaum.github.io/bosc/reference/group_tfr_test.md)
  : Grouped time-frequency map with resampled null testing
- [`group_phase_consistency()`](https://bbuchsbaum.github.io/bosc/reference/group_phase_consistency.md)
  : Grouped phase consistency over bins
- [`group_phase_difference()`](https://bbuchsbaum.github.io/bosc/reference/group_phase_difference.md)
  : Grouped phase difference over bins

## Plotting

- [`plot_spectrum()`](https://bbuchsbaum.github.io/bosc/reference/plot_spectrum.md)
  : Plot spectrum
- [`plot_phase_hist()`](https://bbuchsbaum.github.io/bosc/reference/plot_phase_hist.md)
  : Plot phase histogram
- [`plot_simplebox()`](https://bbuchsbaum.github.io/bosc/reference/plot_simplebox.md)
  : Simple Box Plot by Group
- [`plot_cluster_heatmap()`](https://bbuchsbaum.github.io/bosc/reference/plot_cluster_heatmap.md)
  : Plot Cluster Heatmap

## I/O

- [`read_mat()`](https://bbuchsbaum.github.io/bosc/reference/read_mat.md)
  : Load MATLAB .mat files (wrapper)
- [`open_hdf5()`](https://bbuchsbaum.github.io/bosc/reference/open_hdf5.md)
  : Load HDF5 files (wrapper)
