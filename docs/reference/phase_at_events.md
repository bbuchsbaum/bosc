# Extract narrowband phase at event times

Builds a continuous response trace from event times, narrowband filters
it, computes the analytic signal, and returns the instantaneous phase at
each event. Optionally computes leave-one-out phases (excluding each
event from the trace).

## Usage

``` r
phase_at_events(
  events,
  dt,
  freqlim = NULL,
  fosc = NULL,
  bandwidth = 0.5,
  sd_smooth = NULL,
  leave_one_out = FALSE,
  hilbert_tol = 100,
  filtorder = 2,
  demean = TRUE,
  ...
)
```

## Arguments

- events:

  numeric vector of event times (seconds).

- dt:

  time step for continuous trace (seconds).

- freqlim:

  length-2 numeric vector of filter band. If NULL, supply `fosc`.

- fosc:

  optional center frequency (Hz) used to derive `freqlim`.

- bandwidth:

  half-width (Hz) around `fosc` if `freqlim` is NULL.

- sd_smooth:

  smoothing SD for trace (seconds). If NULL and `fosc` provided,
  defaults to `1/(fosc*8)`, a narrowband smoothing heuristic.

- leave_one_out:

  logical; if TRUE compute phases per event from traces excluding that
  event.

- hilbert_tol, filtorder, demean:

  passed to [`narrowband_hilbert`](narrowband_hilbert.md).

- ...:

  further arguments passed to
  [`make_continuous_trace`](make_continuous_trace.md).

## Value

numeric vector of phases (radians), same length as `events`.
