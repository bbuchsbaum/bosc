#' Narrowband filtering and Hilbert transform
#'
#' Filters data in a narrow frequency band and computes the analytic signal via
#' the Hilbert transform. Mirrors MATLAB \code{narrowbandHilbert} with stability
#' checks and optional demeaning.
#'
#' @param data numeric vector.
#' @param fs sampling rate (Hz).
#' @param freqlim length-2 numeric vector \code{c(fmin, fmax)} in Hz.
#' @param tol tolerance for filter stability (max absolute value before flagging).
#' @param filtorder filter order (integer).
#' @param demean logical; remove mean before filtering.
#' @return list with \code{analytic} (complex), \code{filtered} (numeric).
#' @export
narrowband_hilbert <- function(data,
                               fs,
                               freqlim,
                               tol = 1e2,
                               filtorder = 2,
                               demean = TRUE) {
  if (length(fs) != 1 || !is.finite(fs) || fs <= 0) stop("fs must be positive scalar.")
  if (length(freqlim) != 2) stop("freqlim must be length 2.")
  if (length(filtorder) != 1 || filtorder <= 0) stop("filtorder must be positive.")
  x <- as.numeric(data)
  if (demean) x <- x - mean(x, na.rm = TRUE)

  if (freqlim[1] <= 0) {
    ba <- signal::butter(filtorder * 2, (2 * freqlim[2]) / fs, type = "low")
  } else {
    ba <- signal::butter(filtorder, 2 * freqlim / fs, type = "pass")
  }

  xf <- as.numeric(signal::filtfilt(ba, x))
  if (any(abs(xf) > tol, na.rm = TRUE)) {
    warning("Filter is unstable: values exceed tolerance.")
    return(list(analytic = rep(NA_real_, length(x)), filtered = xf))
  }

  analytic <- analytic_signal(xf - mean(xf, na.rm = TRUE))
  list(analytic = analytic, filtered = xf)
}

# minimal analytic signal via FFT (no external dependency)
analytic_signal <- function(x) {
  n <- length(x)
  X <- fft(x)
  h <- numeric(n)
  if (n %% 2 == 0) {
    h[1] <- 1
    h[n / 2 + 1] <- 1
    h[2:(n / 2)] <- 2
  } else {
    h[1] <- 1
    h[2:((n + 1) / 2)] <- 2
  }
  fft(X * h, inverse = TRUE) / n
}
