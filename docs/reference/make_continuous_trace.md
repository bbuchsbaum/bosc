# Make a continuous trace from event times

Converts discrete event times to a continuous time series with optional
smoothing and quantile trimming (behavior inspired by the MATLAB
`makeContinuousTrace`).

## Usage

``` r
make_continuous_trace(
  events,
  dt,
  sd_smooth = NULL,
  width_block = NULL,
  remove_val = NULL,
  quantlim = NULL,
  warn = TRUE
)
```

## Arguments

- events:

  numeric vector of event times (seconds).

- dt:

  time step for the output series (seconds, \> 0).

- sd_smooth:

  optional standard deviation of a Gaussian kernel (seconds).

- width_block:

  optional width of a block (boxcar) kernel (seconds).

- remove_val:

  optional value to drop (e.g., trial markers) within `dt`.

- quantlim:

  optional length-2 numeric vector of quantile bounds to trim the trace.

- warn:

  logical; emit warnings on dropped points.

## Value

list with `signal` (numeric vector) and `tspan` (time axis).
