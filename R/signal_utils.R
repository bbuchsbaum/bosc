#' @importFrom stats fft
#' @importFrom utils tail
NULL

#' Make a continuous trace from event times
#'
#' Converts discrete event times to a continuous time series with optional smoothing
#' and quantile trimming (behavior inspired by the MATLAB `makeContinuousTrace`).
#'
#' @param events numeric vector of event times (seconds).
#' @param dt time step for the output series (seconds, > 0).
#' @param sd_smooth optional standard deviation of a Gaussian kernel (seconds).
#' @param width_block optional width of a block (boxcar) kernel (seconds).
#' @param remove_val optional value to drop (e.g., trial markers) within \code{dt}.
#' @param quantlim optional length-2 numeric vector of quantile bounds to trim the trace.
#' @param warn logical; emit warnings on dropped points.
#' @return list with \code{signal} (numeric vector) and \code{tspan} (time axis).
#' @export
make_continuous_trace <- function(events,
                                  dt,
                                  sd_smooth = NULL,
                                  width_block = NULL,
                                  remove_val = NULL,
                                  quantlim = NULL,
                                  warn = TRUE) {
  stopifnot(is.numeric(events))
  if (length(dt) != 1 || !is.finite(dt) || dt <= 0) {
    stop("dt must be a positive scalar.")
  }
  if (length(events) == 0) {
    return(list(signal = numeric(0), tspan = numeric(0)))
  }

  ev <- as.numeric(events)
  if (!is.null(remove_val)) {
    keep <- abs(ev - remove_val) > dt
    if (warn && any(!keep)) warning("Dropped ", sum(!keep), " events at remove_val.")
    ev <- ev[keep]
  }
  ev <- ev[!is.na(ev)]
  if (length(ev) == 0) {
    return(list(signal = numeric(0), tspan = numeric(0)))
  }

  tspan <- seq(0, max(ev) + dt, by = dt)
  counts <- numeric(length(tspan))
  idx <- pmax(1L, pmin(length(tspan), round(ev / dt) + 1L))
  for (i in idx) counts[i] <- counts[i] + 1

  if (!is.null(quantlim)) {
    if (length(quantlim) != 2) stop("quantlim must be length 2.")
    nz <- which(counts > 0)
    if (length(nz) > 0) {
      qs <- as.integer(stats::quantile(nz, probs = quantlim, type = 1))
      qs[1] <- max(1L, qs[1])
      qs[2] <- min(length(counts), qs[2])
      counts <- counts[qs[1]:qs[2]]
      tspan <- tspan[qs[1]:qs[2]]
    }
  }

  kernel <- build_kernel(dt, sd_smooth = sd_smooth, width_block = width_block)
  signal <- conv_same(counts, kernel)

  if (warn && sum(counts) != length(ev)) {
    warning("Some events were not mapped to output bins.")
  }

  list(signal = signal, tspan = tspan)
}

#' Centered autocorrelation (unnormalized)
#'
#' Computes an unnormalized autocorrelation with zero lag at the center,
#' equivalent to MATLAB \code{xcorr(x)}.
#'
#' @param x numeric vector.
#' @return numeric vector of length \code{2*length(x) - 1}.
#' @export
autocorr_centered <- function(x) {
  n <- length(x)
  if (n == 0) return(numeric(0))
  lags <- seq(-(n - 1), n - 1)
  res <- numeric(length(lags))
  for (i in seq_along(lags)) {
    lag <- lags[i]
    if (lag >= 0) {
      res[i] <- sum(x[(1 + lag):n] * x[1:(n - lag)])
    } else {
      res[i] <- sum(x[1:(n + lag)] * x[(1 - lag):n])
    }
  }
  res
}

#' Find spectral peak within bounds
#'
#' Computes the FFT power spectrum (single-sided), optional 1/f correction,
#' and returns the frequency of the highest peak within bounds.
#'
#' @param data numeric vector.
#' @param fs sampling rate (Hz).
#' @param flim optional length-2 frequency bounds \code{c(fmin, fmax)} in Hz.
#' @param fcor logical; apply crude 1/f correction.
#' @param taper taper to apply; one of \code{"none"}, \code{"hann"}, \code{"hanning"}.
#' @return list with \code{freq} (peak frequency or NA), \code{fxx} (freq axis),
#'   and \code{spectrum} (power values aligned with \code{fxx}).
#' @export
spectral_peak <- function(data,
                          fs,
                          flim = NULL,
                          fcor = FALSE,
                          taper = c("none", "hann", "hanning")) {
  if (length(fs) != 1 || !is.finite(fs) || fs <= 0) stop("fs must be a positive scalar.")
  taper <- match.arg(taper)
  x <- as.numeric(data)
  L <- length(x)
  if (L == 0) return(list(freq = NA_real_, fxx = numeric(0), spectrum = numeric(0)))

  if (taper != "none") {
    tap <- hanning_window(L)
    x <- x * tap
  }

  tmp <- fft(x)
  ps <- Mod(tmp[seq_len(ceiling(L / 2))]) / L
  if (length(ps) > 2) {
    ps[2:(length(ps) - 1)] <- 2 * ps[2:(length(ps) - 1)]
  }
  fx <- seq(0, fs / 2, length.out = length(ps))

  if (!is.null(flim)) {
    if (length(flim) != 2) stop("flim must be length 2.")
    id <- which(fx >= flim[1] & fx <= flim[2])
    fxx <- fx[id]
    datfft <- ps[id]
  } else {
    fxx <- fx
    datfft <- ps
  }

  corfft <- datfft
  if (fcor) {
    valid <- corfft > 0 & fxx > 0
    if (sum(valid) >= 2) {
      fit <- stats::lm(log10(corfft[valid]) ~ log10(fxx[valid]))
      trend <- 10 ^ stats::predict(fit, newdata = data.frame(fxx = fxx))
      corfft <- corfft - trend
    } else {
      warning("Insufficient positive values for 1/f correction; returning uncorrected spectrum.")
    }
  }

  peak_freq <- NA_real_
  if (length(corfft) >= 3) {
    peaks <- pracma::findpeaks(corfft)
    if (!is.null(peaks) && nrow(peaks) > 0) {
      # pracma::findpeaks returns columns: pks, locs, left_ips, right_ips
      loc <- peaks[which.max(peaks[, 1]), 2]
      peak_freq <- fxx[loc]
    }
  }

  list(freq = peak_freq, fxx = fxx, spectrum = corfft)
}

# ---- internal helpers ----

build_kernel <- function(dt, sd_smooth = NULL, width_block = NULL) {
  if (!is.null(sd_smooth)) {
    sd_pts <- round(sd_smooth / dt)
    if (is.na(sd_pts) || sd_pts <= 0) {
      return(c(0, 1, 0))
    }
    gt <- seq_len(8 * sd_pts)
    return(stats::dnorm(gt, mean = 4 * sd_pts, sd = sd_pts))
  }
  if (!is.null(width_block)) {
    width_pts <- max(1L, round(width_block / dt))
    return(rep(1, width_pts))
  }
  c(0, 1, 0)
}

conv_same <- function(x, kernel) {
  if (length(x) == 0 || length(kernel) == 0) return(numeric(0))
  y_full <- stats::convolve(x, rev(kernel), type = "open")
  k <- length(kernel)
  start <- floor((k - 1) / 2) + 1
  end <- start + length(x) - 1
  as.numeric(y_full[start:end])
}

hanning_window <- function(L) {
  if (L <= 1) return(rep(1, L))
  # Symmetric Hann window matching MATLAB's hann(L) / hanning(L)
  # Formula: 0.5 * (1 - cos(2*pi*n/(L-1))) for n = 0,1,...,L-1
  n <- seq(0, L - 1)
  0.5 * (1 - cos(2 * pi * n / (L - 1)))
}
