#' Configuration-list wrapper for oscillation scores
#'
#' Accepts a list of parameters (named like \code{\link{oscillation_score}}) and
#' forwards them to the core O-score function. Useful when your analysis builds
#' configs programmatically.
#'
#' Recognized fields include \code{fs}, \code{flim}, \code{quantlim},
#' \code{smoothach}, \code{smoothwind}, \code{peakwind}, \code{thresangle},
#' \code{mincycles}, \code{minfreqbandwidth}, \code{fpeak}, \code{warnings},
#' \code{fcor}, \code{taper}, and \code{plot}.
#'
#' @param config list of settings.
#' @param signal numeric vector (binary or continuous).
#' @return list as returned by \code{\link{oscillation_score}}.
#' @export
oscillation_score_config <- function(config, signal) {
  if (is.null(config) || !is.list(config)) stop("config must be a list.")
  args <- config

  if (!is.null(args$warnings) && is.character(args$warnings)) {
    args$warnings <- tolower(args$warnings) %in% c("on", "true", "yes", "1")
  }
  if (!is.null(args$taper) && is.character(args$taper)) {
    args$taper <- tolower(args$taper)
  }

  allowed <- c(
    "fs", "flim", "quantlim", "smoothach", "smoothwind", "peakwind",
    "thresangle", "mincycles", "minfreqbandwidth", "fpeak", "warnings",
    "fcor", "taper", "plot"
  )
  args <- args[names(args) %in% allowed]

  do.call(oscillation_score, c(list(signal = signal), args))
}

#' Configuration-list wrapper for surrogate O-scores
#'
#' Companion to \code{\link{oscillation_score_config}} for surrogate testing.
#' Accepts a config list (named like \code{\link{oscillation_score_surrogates}})
#' and forwards it to the surrogate function.
#'
#' Recognized fields include \code{fs}, \code{flim}, \code{nrep}, \code{fpeak},
#' \code{keep_trend}, \code{trend_dist}, \code{trend_ddt}, \code{trend_alpha},
#' \code{warnings}, \code{fcor}, and \code{taper}.
#'
#' @param config list of settings.
#' @param signal numeric vector.
#' @return list as returned by \code{\link{oscillation_score_surrogates}}.
#' @export
oscillation_score_surrogates_config <- function(config, signal) {
  if (is.null(config) || !is.list(config)) stop("config must be a list.")
  args <- config

  if (!is.null(args$warnings) && is.character(args$warnings)) {
    args$warnings <- tolower(args$warnings) %in% c("on", "true", "yes", "1")
  }
  if (!is.null(args$trend_dist)) {
    args$trend_dist <- tolower(args$trend_dist)
  }
  if (!is.null(args$taper) && is.character(args$taper)) {
    args$taper <- tolower(args$taper)
  }

  allowed <- c(
    "fs", "flim", "nrep", "fpeak", "keep_trend",
    "trend_dist", "trend_ddt", "trend_alpha", "warnings",
    "fcor", "taper"
  )
  args <- args[names(args) %in% allowed]

  do.call(oscillation_score_surrogates, c(list(signal = signal), args))
}

#' @export
#' @rdname oscillation_score_config
oscillation_score_cfg <- function(cfg, signal) {
  .Deprecated("oscillation_score_config")
  oscillation_score_config(cfg, signal)
}

#' @export
#' @rdname oscillation_score_surrogates_config
oscillation_score_stats <- function(cfg, signal) {
  .Deprecated("oscillation_score_surrogates_config")
  oscillation_score_surrogates_config(cfg, signal)
}

#' Oscillation score with surrogate log-Z
#'
#' Computes the oscillation score, generates a surrogate distribution, and returns
#' the log-Z score:
#' \deqn{Z = (log(O) - mean(log(O_{rep}))) / sd(log(O_{rep}))}.
#'
#' @param signal numeric vector (binary or continuous).
#' @param fs sampling rate in Hz.
#' @param flim length-2 numeric vector \code{c(fmin, fmax)} in Hz.
#' @param nrep number of surrogates.
#' @param alpha significance level for one-tailed Z-thresholding.
#' @param keep_trend,trend_dist,trend_ddt,trend_alpha passed to surrogates.
#' @param fcor logical; apply 1/f correction.
#' @param taper taper before FFT.
#' @param ... further arguments passed to \code{\link{oscillation_score}}.
#' @param tidy logical; if TRUE return a one-row data.frame.
#' @return list (default) or data.frame if \code{tidy=TRUE}.
#' @export
oscillation_score_z <- function(signal,
                     fs,
                     flim,
                     nrep = 500,
                     alpha = 0.05,
                     keep_trend = FALSE,
                     trend_dist = c("gamma"),
                     trend_ddt = NULL,
                     trend_alpha = 0.05,
                     fcor = FALSE,
                     taper = c("none", "hann", "hanning"),
                     tidy = FALSE,
                     ...) {
  taper <- match.arg(taper)

  os <- oscillation_score(
    signal = signal,
    fs = fs,
    flim = flim,
    warnings = FALSE,
    fcor = fcor,
    taper = taper,
    ...
  )

  if (!is.finite(os$oscore) || os$oscore <= 0 || !is.finite(os$fosc)) {
    out <- list(
      oscore = os$oscore,
      fosc = os$fosc,
      surrogates = NULL,
      z = NA_real_,
      pval = NA_real_,
      significant = FALSE,
      flim = os$flim
    )
    if (isTRUE(tidy)) return(as.data.frame(out[setdiff(names(out), "surrogates")]))
    return(out)
  }

  sur <- oscillation_score_surrogates(
    signal = signal,
    fs = fs,
    flim = flim,
    nrep = nrep,
    fpeak = os$fosc,
    keep_trend = keep_trend,
    trend_dist = trend_dist,
    trend_ddt = trend_ddt,
    trend_alpha = trend_alpha,
    warnings = FALSE,
    fcor = fcor,
    taper = taper
  )

  valid <- sur$oscore_rp
  valid <- valid[is.finite(valid) & valid > 0]

  z <- NA_real_
  pval <- NA_real_
  significant <- FALSE
  if (length(valid) >= 2) {
    z <- (log(os$oscore) - mean(log(valid))) / stats::sd(log(valid))
    pval <- 1 - stats::pnorm(z)
    significant <- is.finite(z) && z >= stats::qnorm(1 - alpha)
  }

  out <- list(
    oscore = os$oscore,
    fosc = os$fosc,
    surrogates = sur,
    z = z,
    pval = pval,
    significant = significant,
    flim = os$flim
  )
  if (isTRUE(tidy)) {
    return(as.data.frame(out[setdiff(names(out), "surrogates")]))
  }
  out
}

#' @export
#' @rdname oscillation_score_z
oscore_z <- function(...) {
  .Deprecated("oscillation_score_z")
  oscillation_score_z(...)
}

#' Extract narrowband phase at event times
#'
#' Builds a continuous response trace from event times, narrowband filters it,
#' computes the analytic signal, and returns the instantaneous phase at each event.
#' Optionally computes leave-one-out phases (excluding each event from the trace).
#'
#' @param events numeric vector of event times (seconds).
#' @param dt time step for continuous trace (seconds).
#' @param freqlim length-2 numeric vector of filter band. If NULL, supply \code{fosc}.
#' @param fosc optional center frequency (Hz) used to derive \code{freqlim}.
#' @param bandwidth half-width (Hz) around \code{fosc} if \code{freqlim} is NULL.
#' @param sd_smooth smoothing SD for trace (seconds). If NULL and \code{fosc} provided,
#'   defaults to \code{1/(fosc*8)}, a narrowband smoothing heuristic.
#' @param leave_one_out logical; if TRUE compute phases per event from traces
#'   excluding that event.
#' @param hilbert_tol,filtorder,demean passed to \code{\link{narrowband_hilbert}}.
#' @param ... further arguments passed to \code{\link{make_continuous_trace}}.
#' @return numeric vector of phases (radians), same length as \code{events}.
#' @export
phase_at_events <- function(events,
                            dt,
                            freqlim = NULL,
                            fosc = NULL,
                            bandwidth = 0.5,
                            sd_smooth = NULL,
                            leave_one_out = FALSE,
                            hilbert_tol = 1e2,
                            filtorder = 2,
                            demean = TRUE,
                            ...) {
  if (!is.numeric(events)) stop("events must be numeric.")
  if (length(dt) != 1 || !is.finite(dt) || dt <= 0) stop("dt must be positive scalar.")

  ev_all <- as.numeric(events)
  ev <- ev_all[!is.na(ev_all)]
  if (length(ev) == 0) return(rep(NA_real_, length(events)))

  fs <- 1 / dt

  if (is.null(freqlim)) {
    if (is.null(fosc) || !is.finite(fosc)) stop("Provide freqlim or a finite fosc.")
    freqlim <- c(fosc - bandwidth, fosc + bandwidth)
  }

  if (is.null(sd_smooth) && !is.null(fosc) && is.finite(fosc) && fosc > 0) {
    sd_smooth <- 1 / (fosc * 8)
  }

  compute_phases <- function(ref_events) {
    trace <- make_continuous_trace(ref_events, dt = dt, sd_smooth = sd_smooth, warn = FALSE, ...)
    if (length(trace$signal) == 0) return(rep(NA_real_, length(events)))

    nb <- narrowband_hilbert(
      data = trace$signal,
      fs = fs,
      freqlim = freqlim,
      tol = hilbert_tol,
      filtorder = filtorder,
      demean = demean
    )
    phases_all <- Arg(nb$analytic)

    idx <- match(round(ev_all / dt), round(trace$tspan / dt))
    out <- rep(NA_real_, length(events))
    ok <- !is.na(idx) & idx > 0 & idx <= length(phases_all)
    out[ok] <- phases_all[idx[ok]]
    out
  }

  if (!leave_one_out) {
    return(compute_phases(ev))
  }

  out <- rep(NA_real_, length(events))
  for (i in seq_along(events)) {
    if (is.na(events[i])) next
    ref <- ev
    pos <- match(events[i], ref)
    if (!is.na(pos)) ref <- ref[-pos]
    if (length(ref) < 2) next
    out[i] <- compute_phases(ref)[i]
  }
  out
}
