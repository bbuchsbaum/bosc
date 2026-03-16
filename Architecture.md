# Generalized Group Analysis API for `bosc`

## Goal

Incorporate the reusable methodological ideas from `MBO_shared` without
importing paper-specific concepts, data assumptions, or pipeline structure.

The result should feel like a natural extension of the current `bosc` API:

- compact surface area,
- composable low-level functions,
- config/wrapper ergonomics where appropriate,
- stable list return objects with tidy accessors,
- no manuscript-specific vocabulary in exported names.

## Design Principles

1. Port abstractions, not scripts.
2. Keep preprocessing separate from analysis.
3. Make grouped methods work for any binned behavioral series, not only SOA data.
4. Prefer a small number of orthogonal verbs over many convenience wrappers.
5. Use one null-model interface across grouped methods.
6. Reuse existing `bosc` internals wherever possible.

## What Belongs in Scope

The reusable content from `MBO_shared` falls into five generic layers:

1. Trial preprocessing
2. Binning and aggregation
3. Grouped spectral summaries
4. Grouped null-model testing
5. Grouped time-frequency and phase summaries

The following do not belong in the package API:

- nicotine analyses,
- "in the zone" analyses,
- manuscript plotting helpers,
- experiment-specific column names,
- paper-specific contrasts and labels.

## Fit With Current `bosc`

The current package already has strong low-level building blocks:

- signal construction and spectrum helpers,
- oscillation score and surrogate workflows,
- phase extraction helpers,
- circular statistics,
- cluster detection,
- tidy summary accessors.

What is missing is an elegant layer for grouped behavioral analyses over
trial-level data or already-binned subject-level series.

## Proposed Abstraction Layer

### 1. Preprocessing

These functions operate on trial tables and return trial tables or binned tables.

```r
align_time_bins(data, time, bins, by = NULL, tolerance = NULL, method = "nearest")
residualize_trials(data, response, terms = NULL, by = NULL, family = "gaussian")
aggregate_by_bin(data, value, bin, by = NULL, fun = mean, na_rm = TRUE, complete = TRUE)
```

Intent:

- `align_time_bins()` maps observed timing values onto a target grid.
- `residualize_trials()` removes nuisance structure without hard-coding any paper covariates.
- `aggregate_by_bin()` converts trial tables into regularized per-subject series.

These should accept bare vectors or data-frame columns by standard R conventions,
but the first implementation can stay simple and explicit if needed.

### 2. Grouped Series Object

Introduce a lightweight internal convention for binned subject-level data:

```r
list(
  data = matrix_or_array,
  subjects = character(),
  bins = numeric(),
  measure = "hr",
  meta = list(...)
)
```

This does not need a formal S3 class in the first pass, but the structure should
be consistent enough that later methods can dispatch on it cleanly.

### 3. Grouped Spectral API

These functions provide the generic equivalent of the paper's grouped FFT logic.

```r
group_spectrum(x, fs, flim = NULL, detrend = TRUE, taper = "hann", average = TRUE)
group_spectrum_test(
  x,
  fs,
  flim = NULL,
  null = c("shuffle_labels", "sign_flip", "circular_shift", "ar1"),
  nrep = 1000,
  detrend = TRUE,
  taper = "hann",
  p_adjust = c("none", "fdr")
)
```

`x` should accept either:

- a grouped-series object, or
- a numeric matrix with shape `subjects x bins`.

Return shape:

```r
list(
  observed = list(freq = ..., power = ..., per_subject = ...),
  null = list(power = ..., method = ..., nrep = ...),
  stats = list(p = ..., z = ..., p_adj = ..., significant = ...),
  meta = list(fs = ..., flim = ..., detrend = ..., taper = ...)
)
```

### 4. Grouped Time-Frequency API

This generalizes the paper's wavelet/TFA analyses without carrying over its
manuscript-specific plotting layer.

```r
group_tfr(
  x,
  fs,
  freqs,
  method = c("morlet"),
  n_cycles = NULL,
  reflect = FALSE,
  edge = 0
)

group_tfr_test(
  x,
  fs,
  freqs,
  null = c("shuffle_labels", "sign_flip", "circular_shift", "ar1"),
  nrep = 1000,
  reflect = FALSE,
  edge = 0,
  cluster = TRUE
)
```

Return shape:

```r
list(
  observed = list(map = ..., per_subject = ..., time = ..., freq = ...),
  null = list(map = ..., method = ..., nrep = ...),
  stats = list(z = ..., p = ..., clusters = ...),
  meta = list(...)
)
```

Cluster correction should reuse `detect_clusters()` rather than reimplementing
image-label logic from the paper code.

### 5. Grouped Phase API

These functions generalize participant-level phase concentration and
phase-difference summaries.

```r
group_phase_consistency(x, fs, freq, detrend = TRUE, taper = "hann")
group_phase_difference(x1, x2, fs, freq, mode = c("signed_pi", "wrapped_2pi"))
```

Return shape:

```r
list(
  observed = list(phase = ..., summary = ...),
  stats = list(r = ..., p = ..., rayleigh = ...),
  meta = list(freq = ..., fs = ...)
)
```

These should align with the semantics of `phase_at_events()` and current
circular-stat helpers.

## Null Model Design

The grouped-analysis layer should use a shared null-model interface rather than
embedding separate permutation logic in each function.

Proposed internal helpers:

```r
generate_group_null(x, method, ...)
apply_group_null(x, null_draw, statistic, ...)
```

Supported methods, in rollout order:

1. `sign_flip`
2. `shuffle_labels`
3. `circular_shift`
4. `ar1`

Notes:

- `sign_flip` is appropriate for paired/difference spectra.
- `shuffle_labels` is appropriate when labels are exchangeable within strata.
- `circular_shift` is appropriate for preserving within-series autocorrelation.
- `ar1` is valuable, but should be added only after the generic null interface is stable.

The `AR1_control` script in `MBO_shared` should be treated as a conceptual
reference only, not as source code to port directly.

## Naming and Ergonomics

Naming should follow the package's current style:

- verb-first exported function names,
- short argument lists,
- explicit defaults,
- plain list returns,
- optional tidy helpers.

Avoid names like:

- `soa_avg()`
- `FFT_grouped()`
- `fullStatTFA()`
- `grpPhaseCon()`

Prefer:

- `aggregate_by_bin()`
- `group_spectrum()`
- `group_tfr_test()`
- `group_phase_consistency()`

## Tidy Accessors

If the grouped APIs land, the package should also expose matching accessors:

```r
group_spectrum_tidy(x)
group_tfr_tidy(x)
group_phase_tidy(x)
```

These should mirror the existing `oscore_tidy()`, `oscore_spectrum()`,
`clusters_tidy()`, and `clusters_pixels()` conventions.

## Staged Implementation Plan

### Phase 1: Foundation

Implement:

- `align_time_bins()`
- `residualize_trials()`
- `aggregate_by_bin()`
- internal grouped-series validator/helper

Deliverable:

- clean preprocessing path from trial table to grouped matrix.

### Phase 2: Grouped Spectrum

Implement:

- `group_spectrum()`
- `group_spectrum_test()`
- optional `group_spectrum_tidy()`

Deliverable:

- general grouped FFT + null testing workflow.

### Phase 3: Grouped TFR and Phase

Implement:

- `group_tfr()`
- `group_tfr_test()`
- `group_phase_consistency()`
- `group_phase_difference()`

Deliverable:

- coherent wavelet/phase layer for group analyses.

### Phase 4: AR1 Nulls

Implement:

- `ar1` backend for grouped null generation

Deliverable:

- robust autocorrelation-aware null model with tests and clear constraints.

## Test Strategy

All new computational features should ship with explicit contracts and at least
four test families:

1. Contract tests
2. Metamorphic tests
3. Differential/oracle tests
4. Edge/adversarial tests

### Phase 1 tests

- bin alignment is stable under exact-grid and near-grid inputs,
- residualization leaves nuisance-free synthetic data unchanged,
- aggregation preserves subject/bin completeness rules,
- missing values and empty strata fail predictably.

### Phase 2 tests

- grouped spectra recover known oscillations in simulated subject matrices,
- null p-values are approximately uniform under null simulations,
- sign-flip and shuffle paths respect shape and finite-value invariants,
- grouped results match direct per-subject calls where mathematically equivalent.

### Phase 3 tests

- TFR peak localization tracks known injected time-frequency structure,
- cluster detection is invariant to equivalent matrix representations,
- phase consistency increases when subject phases are aligned,
- phase-difference summaries behave correctly under known phase offsets.

### Phase 4 tests

- AR1 nulls preserve approximate lag-1 structure,
- false-positive rates are controlled on autocorrelated null data,
- AR1 path degrades gracefully on short or degenerate series.

## Recommended Immediate Next Step

Start with Phase 1 and Phase 2 only.

That gives `bosc` a general behavioral group-analysis layer quickly, without
locking in premature decisions about wavelet internals or AR1 implementation
details.
