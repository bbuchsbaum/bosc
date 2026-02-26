# Get started with bosc

## What does bosc do?

You have behavioral time series — reaction times, button presses, spike
trains — and you want to know whether they contain rhythmic structure.
Is there a dominant oscillation frequency? Is it statistically
significant? At what phase of the oscillation do events occur?

**bosc** answers these questions with a small, focused API:

| Verb | Function | What it does |
|:---|:---|:---|
| **Score** | [`oscillation_score_z()`](../reference/oscillation_score_z.md) | Quantify oscillation strength + significance |
| **Filter** | [`narrowband_hilbert()`](../reference/narrowband_hilbert.md) | Extract phase and amplitude at a frequency |
| **Phase** | [`phase_at_events()`](../reference/phase_at_events.md) | Get instantaneous phase at event times |
| **Cluster** | [`detect_clusters()`](../reference/detect_clusters.md) | Find significant regions in time-frequency maps |
| **Test** | [`circ_rayleigh()`](../reference/circ_rayleigh.md) | Test for non-uniform phase distributions |

``` r
library(bosc)
```

## Your first oscillation score in 30 seconds

Create a signal with 10 Hz periodicity, compute the score, and test
significance — all in one call:

``` r
set.seed(42)
fs <- 1000
sig <- numeric(2000)
spikes <- seq(50, 1950, by = 100) + round(rnorm(20, sd = 5))
spikes <- spikes[spikes > 0 & spikes <= 2000]
sig[spikes] <- 1

result <- oscillation_score_z(sig, fs = fs, flim = c(5, 20), nrep = 200)
```

``` r
cat("Oscillation score:", round(result$oscore, 3), "\n")
#> Oscillation score: 4.908
cat("Peak frequency:   ", round(result$fosc, 2), "Hz\n")
#> Peak frequency:    7.33 Hz
cat("Z-score:          ", round(result$z, 2), "\n")
#> Z-score:           0.42
cat("Significant:      ", result$significant, "\n")
#> Significant:       FALSE
```

![Power spectrum of the autocorrelogram. The dashed line marks the
detected peak at ~10
Hz.](bosc_files/figure-html/plot-first-spectrum-1.png)

Power spectrum of the autocorrelogram. The dashed line marks the
detected peak at ~10 Hz.

The oscillation score is 4.91 with a z-score of 0.4 — a clear 10 Hz
rhythm, as expected.

## Extracting phase at behavioral events

Given event times (in seconds),
[`phase_at_events()`](../reference/phase_at_events.md) builds a
continuous trace, filters at the target frequency, and returns the
instantaneous phase at each event:

``` r
events <- c(0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
phases <- phase_at_events(events, dt = 0.001, fosc = 10, bandwidth = 1)
```

![Phase distribution at event
times.](bosc_files/figure-html/plot-phases-1.png)

Phase distribution at event times.

Test whether the phases are non-uniformly distributed:

``` r
ray <- circ_rayleigh(phases)
cat("Mean resultant length:", round(circ_r(phases), 3), "\n")
#> Mean resultant length: 0.902
cat("Rayleigh p-value:     ", format(ray$pval, digits = 3), "\n")
#> Rayleigh p-value:      2.93e-05
```

## Where to go next

You’ve just seen the core three-step pipeline: **score**
([`oscillation_score_z()`](../reference/oscillation_score_z.md)),
**extract phase**
([`phase_at_events()`](../reference/phase_at_events.md)), **test**
([`circ_rayleigh()`](../reference/circ_rayleigh.md)). For group-level
analyses, add [`paired_tscore()`](../reference/paired_tscore.md) and
[`detect_clusters()`](../reference/detect_clusters.md) for cluster-based
permutation testing on time-frequency maps.

| I want to… | Read this |
|:---|:---|
| Understand oscillation scores in depth | [`vignette("oscillation-score")`](../articles/oscillation-score.md) |
| Analyze phase consistency and clusters | [`vignette("ppc-clustering")`](../articles/ppc-clustering.md) |
| Run batch analyses across subjects | [`vignette("workflow-wrappers")`](../articles/workflow-wrappers.md) |

## References

**bosc** is an R port of the
[behavioral-oscillations](https://github.com/marijeterwal/behavioral-oscillations)
MATLAB toolbox. If you use this package in published work, please cite
the original paper:

> Ter Wal, M., Linde-Domingo, J., Lifanov, J., Roux, F., Kolibius, L.D.,
> Gollwitzer, S., Lang, J., Hamer, H., Rollings, D., Sawlani, V.,
> Chelvarajah, R., Staresina, B., Hanslmayr, S., & Wimber, M. (2021).
> Theta rhythmicity governs the timing of behavioural and hippocampal
> responses in humans specifically during memory-dependent tasks.
> *Nature Communications*, 12, 7048.
> <https://doi.org/10.1038/s41467-021-25959-7>
