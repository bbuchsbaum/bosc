# Package index

## Signal and O-score

- [`make_continuous_trace()`](make_continuous_trace.md) : Make a
  continuous trace from event times
- [`autocorr_centered()`](autocorr_centered.md) : Centered
  autocorrelation (unnormalized)
- [`spectral_peak()`](spectral_peak.md) : Find spectral peak within
  bounds
- [`oscillation_score()`](oscillation_score.md) : Oscillation Score for
  binary or continuous data
- [`oscillation_score_surrogates()`](oscillation_score_surrogates.md) :
  Surrogate oscillation scores (shuffle or trend-preserving)

## Workflow Wrappers

- [`oscillation_score_z()`](oscillation_score_z.md)
  [`oscore_z()`](oscillation_score_z.md) : Oscillation score with
  surrogate log-Z
- [`oscillation_score_config()`](oscillation_score_config.md)
  [`oscillation_score_cfg()`](oscillation_score_config.md) :
  Configuration-list wrapper for oscillation scores
- [`oscillation_score_surrogates_config()`](oscillation_score_surrogates_config.md)
  [`oscillation_score_stats()`](oscillation_score_surrogates_config.md)
  : Configuration-list wrapper for surrogate O-scores
- [`phase_at_events()`](phase_at_events.md) : Extract narrowband phase
  at event times

## Tidy Helpers

- [`oscore_tidy()`](oscore_tidy.md) : Tidy summary of oscillation score
  results
- [`oscore_spectrum()`](oscore_spectrum.md) : Tidy spectrum from
  oscillation score or spectral peak output
- [`clusters_tidy()`](clusters_tidy.md) : Tidy summary of detected
  clusters
- [`clusters_pixels()`](clusters_pixels.md) : Tidy per-pixel cluster
  data

## Phase and Hilbert

- [`narrowband_hilbert()`](narrowband_hilbert.md) : Narrowband filtering
  and Hilbert transform

## Circular Statistics

- [`circ_mean()`](circ_mean.md) : Circular Mean Direction
- [`circ_r()`](circ_r.md) : Mean Resultant Vector Length
- [`circ_rayleigh()`](circ_rayleigh.md) : Rayleigh Test for
  Non-uniformity
- [`circ_vtest()`](circ_vtest.md) : V-test for Non-uniformity with
  Specified Mean Direction

## PPC and Statistical Testing

- [`pairwise_phase_consistency()`](pairwise_phase_consistency.md) :
  Pairwise Phase Consistency (1D)
- [`u_score_matrix()`](u_score_matrix.md) : Mann-Whitney U Z-score per
  element vs. reference repeats
- [`paired_tscore()`](paired_tscore.md) : Paired t-score along a
  dimension
- [`nonparam_pval()`](nonparam_pval.md) : Non-parametric p-values vs
  reference
- [`fdr_bh()`](fdr_bh.md) : FDR Correction (Benjamini-Hochberg)

## Clustering

- [`extract_clusters()`](extract_clusters.md) : Connected component
  labeling for thresholded data
- [`detect_clusters()`](detect_clusters.md) : Detect clusters with
  optional splitting

## Plotting

- [`plot_spectrum()`](plot_spectrum.md) : Plot spectrum
- [`plot_phase_hist()`](plot_phase_hist.md) : Plot phase histogram
- [`plot_simplebox()`](plot_simplebox.md) : Simple Box Plot by Group
- [`plot_cluster_heatmap()`](plot_cluster_heatmap.md) : Plot Cluster
  Heatmap

## I/O

- [`read_mat()`](read_mat.md) : Load MATLAB .mat files (wrapper)
- [`open_hdf5()`](open_hdf5.md) : Load HDF5 files (wrapper)
