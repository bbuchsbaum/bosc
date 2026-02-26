#' Oscillation Score for binary or continuous data
#'
#' Port of the MATLAB \code{oscillationScore} algorithm with robustness fixes.
#' Computes an autocorrelogram, removes the central peak, finds the spectral peak,
#' and returns the oscillation score.
#'
#' @param signal numeric vector (binary or continuous).
#' @param fs sampling rate in Hz.
#' @param flim length-2 numeric vector \code{c(fmin, fmax)} in Hz.
#' @param quantlim optional length-2 quantile bounds for trimming nonzero indices.
#' @param smoothach logical; smooth the ACH with a Gaussian kernel.
#' @param smoothwind smoothing kernel width (seconds); default matches MATLAB heuristic.
#' @param peakwind slow-kernel width for peak removal (seconds); default matches MATLAB heuristic.
#' @param thresangle angular threshold for peak detection (degrees).
#' @param mincycles minimum cycles to include at \code{fmin}.
#' @param minfreqbandwidth optional minimum bandwidth; returns NA if unmet.
#' @param fpeak optional frequency of interest; overrides peak search.
#' @param warnings logical; emit warnings.
#' @param fcor logical; apply 1/f correction in spectral peak estimation.
#' @param taper taper applied before FFT; one of \code{"none"}, \code{"hann"}, \code{"hanning"}.
#' @param signal_mode one of \code{"auto"}, \code{"event"}, or \code{"continuous"}.
#'   \code{"event"} uses event-count heuristics for effective frequency bounds;
#'   \code{"continuous"} uses signal-variance heuristics and Nyquist-limited bounds.
#' @param plot logical; if TRUE, plots intermediate spectra (optional, requires ggplot2).
#' @return list with \code{oscore}, \code{fosc}, \code{flim}, \code{flimfft}, \code{freqs}.
#' @export
oscillation_score <- function(signal,
                              fs,
                              flim,
                              quantlim = NULL,
                              smoothach = TRUE,
                              smoothwind = NULL,
                              peakwind = NULL,
                              thresangle = 10,
                              mincycles = 3,
                              minfreqbandwidth = NULL,
                              fpeak = NULL,
                              warnings = TRUE,
                              fcor = FALSE,
                              taper = c("none", "hann", "hanning"),
                              signal_mode = c("auto", "event", "continuous"),
                              plot = FALSE) {
  stopifnot(is.numeric(signal))
  if (length(fs) != 1 || !is.finite(fs) || fs <= 0) stop("fs must be positive scalar.")
  if (length(flim) != 2) stop("flim must be length 2.")
  if (any(!is.finite(flim))) stop("flim must contain finite values.")
  if (flim[1] <= 0 || flim[2] <= 0 || flim[1] >= flim[2]) {
    stop("flim must be increasing positive frequencies.")
  }
  taper <- match.arg(taper)
  signal_mode <- match.arg(signal_mode)

  sig <- as.numeric(signal)
  sig[is.na(sig)] <- 0

  if (signal_mode == "auto") {
    signal_mode <- if (all(sig >= 0, na.rm = TRUE)) "event" else "continuous"
  }

  if (signal_mode == "event") {
    if (sum(sig, na.rm = TRUE) < 3) {
      if (warnings) warning("Dataset is (nearly) empty: Oscillation score could not be computed.")
      return(list(oscore = NA_real_, fosc = NA_real_, flim = numeric(0), flimfft = numeric(0), freqs = numeric(0)))
    }
  } else {
    if (!is.finite(stats::sd(sig, na.rm = TRUE)) || stats::sd(sig, na.rm = TRUE) < .Machine$double.eps) {
      if (warnings) warning("Signal has near-zero variance: Oscillation score could not be computed.")
      return(list(oscore = NA_real_, fosc = NA_real_, flim = numeric(0), flimfft = numeric(0), freqs = numeric(0)))
    }
  }

  if (!is.null(quantlim)) {
    nz <- which(sig != 0)
    if (length(nz) == 0) {
      if (warnings) warning("No nonzero samples after quantile trimming.")
      return(list(oscore = NA_real_, fosc = NA_real_, flim = numeric(0), flimfft = numeric(0), freqs = numeric(0)))
    }
    q <- as.integer(stats::quantile(nz, probs = quantlim, type = 1))
    q[1] <- max(1L, q[1] - 1L)
    q[2] <- min(length(sig), q[2] + 1L)
    sigq <- sig[q[1]:q[2]]
  } else {
    tstart <- max(1L, which(sig != 0)[1] - 1L)
    tend <- min(length(sig), tail(which(sig != 0), 1) + 1L)
    sigq <- sig[tstart:tend]
  }

  if (length(sigq) == 0) {
    if (warnings) warning("Empty trimmed signal.")
    return(list(oscore = NA_real_, fosc = NA_real_, flim = numeric(0), flimfft = numeric(0), freqs = numeric(0)))
  }

  fmin <- max(flim[1], mincycles * fs / length(sigq))
  if (signal_mode == "event") {
    fmax <- min(flim[2], sum(sigq, na.rm = TRUE) / (length(sigq) / fs))
  } else {
    fmax <- min(flim[2], fs / 2)
  }
  flim_eff <- c(fmin, fmax)

  if (!all(is.finite(flim_eff)) || fmax <= fmin) {
    if (warnings) warning("Effective frequency bounds are invalid for this signal.")
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = numeric(0), freqs = numeric(0)))
  }

  if (!is.null(minfreqbandwidth) && (fmax - fmin) < minfreqbandwidth) {
    if (warnings) warning("Data not sufficient to meet minimal frequency bandwidth.")
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = numeric(0), freqs = numeric(0)))
  }

  w <- 2 ^ (1 + floor(max(log2(2 * mincycles * fs / flim[1]), log2(fs / 2))))

  if (is.null(smoothwind)) smoothwind <- min(0.002, 134 / (1.5 * flim[2]) / 1000)
  if (is.null(peakwind)) peakwind <- 2 * 134 / (1.5 * flim[1]) / 1000

  ach <- autocorr_centered(sigq)

  if (smoothach) {
    sdwind <- round(smoothwind * fs)
    if (sdwind > 0) {
      gt <- seq(-4 * sdwind, 4 * sdwind)
      wind <- (1 / (sdwind * sqrt(2 * pi))) * exp(-(gt ^ 2) / (2 * sdwind ^ 2))
      ach_smooth <- conv_same(ach, wind)
    } else {
      ach_smooth <- ach
    }
  } else {
    ach_smooth <- ach
  }

  sdwind <- round(peakwind * fs)
  if (sdwind > 0) {
    gt <- seq(-4 * sdwind, 4 * sdwind)
    wind <- (1 / (sdwind * sqrt(2 * pi))) * exp(-(gt ^ 2) / (2 * sdwind ^ 2))
    ach_slow <- conv_same(ach, wind)
  } else {
    ach_slow <- ach
  }

  center_idx <- ceiling(length(ach_slow) / 2)
  denom <- ach_slow[center_idx]
  if (is.na(denom) || abs(denom) < .Machine$double.eps) {
    if (warnings) warning("ACH center near zero; cannot determine peak.")
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = numeric(0), freqs = numeric(0)))
  }
  thres <- tan(pi * thresangle / 180)
  scalfact <- (length(ach) - 1) / denom

  half <- ach_slow[seq_len(center_idx)]
  diff_ach <- rev(diff(half))
  if (length(diff_ach) == 0) {
    if (warnings) warning("Dataset is (nearly) empty: Oscillation score could not be computed.")
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = numeric(0), freqs = numeric(0)))
  }

  curvature <- c(diff(diff_ach), 0)
  peak_candidates <- which(scalfact * diff_ach <= thres & curvature < 0)
  phw <- if (length(peak_candidates) == 0) 1 else peak_candidates[1]

  start <- center_idx + phw
  if (start > length(ach_smooth)) {
    if (warnings) warning("Insufficient data after peak removal.")
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = numeric(0), freqs = numeric(0)))
  }
  tmp <- ach_smooth[start:length(ach_smooth)]
  if (length(tmp) < w) {
    ach_nopeak <- numeric(w)
    ach_nopeak[seq_along(tmp)] <- tmp
  } else {
    ach_nopeak <- tmp[seq_len(w)]
  }

  sp <- spectral_peak(ach_nopeak, fs = fs, flim = flim_eff, fcor = fcor, taper = taper)
  peakfreq <- sp$freq
  freqs <- sp$fxx
  flimfft <- sp$spectrum

  if (!is.null(fpeak) && is.finite(fpeak)) {
    id <- which.min(abs(freqs - fpeak))
    peakfreq <- freqs[id]
  }

  tot <- spectral_peak(ach_nopeak, fs = fs, flim = c(0, fs / 2), fcor = fcor, taper = taper)
  if (is.na(peakfreq) || length(freqs) == 0) {
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = flimfft, freqs = freqs))
  }
  id <- which.min(abs(freqs - peakfreq))
  if (length(id) == 0 || is.na(id)) {
    return(list(oscore = NA_real_, fosc = NA_real_, flim = flim_eff, flimfft = flimfft, freqs = freqs))
  }
  denom_mean <- mean(tot$spectrum, na.rm = TRUE)
  if (!is.finite(denom_mean) || denom_mean == 0) {
    if (warnings) warning("Mean spectrum is zero/invalid.")
    return(list(oscore = NA_real_, fosc = peakfreq, flim = flim_eff, flimfft = flimfft, freqs = freqs))
  }
  oscore <- flimfft[id] / denom_mean

  if (isTRUE(plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    try(plot_spectrum(freqs, flimfft, peak = peakfreq), silent = TRUE)
  }

  list(oscore = oscore, fosc = peakfreq, flim = flim_eff, flimfft = flimfft, freqs = freqs)
}

#' Surrogate oscillation scores (shuffle or trend-preserving)
#'
#' Generates surrogate signals and computes oscillation scores for reference distributions.
#' Fixes edge cases from the MATLAB implementation (e.g., negative p-values, missing sums).
#'
#' @param signal numeric vector.
#' @param fs sampling rate in Hz.
#' @param flim length-2 numeric vector \code{c(fmin, fmax)}.
#' @param nrep number of surrogates.
#' @param fpeak peak frequency of original signal (Hz); required for windowing.
#' @param keep_trend logical; if TRUE fit trend distribution instead of shuffling.
#' @param trend_dist character vector of distributions for \code{fitdistrplus::fitdist}.
#' @param trend_ddt step (seconds) for trend generation; defaults to \code{0.5/fs}.
#' @param trend_alpha p-value threshold to accept fit.
#' @param surrogate_method surrogate generator; one of \code{"auto"},
#'   \code{"event_shuffle"}, \code{"event_trend"}, or \code{"phase_randomized"}.
#'   In \code{"auto"} mode, continuous signals default to phase randomization,
#'   while event-like signals use shuffle or trend fitting depending on \code{keep_trend}.
#' @param signal_mode one of \code{"auto"}, \code{"event"}, or \code{"continuous"}.
#' @param warnings logical; emit warnings.
#' @param fcor logical; apply 1/f correction in surrogate peak estimation.
#' @param taper taper applied before FFT; one of \code{"none"}, \code{"hann"}, \code{"hanning"}.
#' @return list with \code{oscore_rp}, \code{fosc_rp}, \code{trendfit},
#'   \code{surrogate_method}, \code{signal_mode}, and (optionally) \code{signrep}.
#' @export
oscillation_score_surrogates <- function(signal,
                                         fs,
                                         flim,
                                         nrep,
                                         fpeak,
                                         keep_trend = FALSE,
                                         trend_dist = c("gamma"),
                                         trend_ddt = NULL,
                                         trend_alpha = 0.05,
                                         surrogate_method = c(
                                           "auto", "event_shuffle", "event_trend", "phase_randomized"
                                         ),
                                         signal_mode = c("auto", "event", "continuous"),
                                         warnings = TRUE,
                                         fcor = FALSE,
                                         taper = c("none", "hann", "hanning")) {
  if (length(nrep) != 1 || !is.finite(nrep) || nrep <= 0) stop("nrep must be positive.")
  if (length(fs) != 1 || !is.finite(fs) || fs <= 0) stop("fs must be positive.")
  if (length(flim) != 2 || any(!is.finite(flim)) || flim[1] <= 0 || flim[2] <= 0 || flim[1] >= flim[2]) {
    stop("flim must be length 2 with increasing positive values.")
  }
  if (is.null(trend_ddt)) trend_ddt <- 0.5 / fs
  taper <- match.arg(taper)
  surrogate_method <- match.arg(surrogate_method)
  signal_mode <- match.arg(signal_mode)

  sig <- as.numeric(signal)
  sig[is.na(sig)] <- 0
  if (signal_mode == "auto") {
    signal_mode <- if (all(sig >= 0, na.rm = TRUE)) "event" else "continuous"
  }

  if (surrogate_method == "auto") {
    surrogate_method <- if (signal_mode == "continuous") {
      "phase_randomized"
    } else if (isTRUE(keep_trend)) {
      "event_trend"
    } else {
      "event_shuffle"
    }
  }
  if (isTRUE(keep_trend) && surrogate_method != "event_trend" && warnings) {
    warning("keep_trend ignored because surrogate_method != 'event_trend'.")
  }

  if (surrogate_method != "phase_randomized" && is.null(fpeak)) {
    stop("fpeak must be provided for event-based surrogate methods (can be NA to use flim[1] window).")
  }

  if (signal_mode == "event") {
    if (sum(sig, na.rm = TRUE) < 3) {
      if (warnings) warning("Oscillation score stats could not be computed.")
      return(list(
        oscore_rp = NA_real_,
        fosc_rp = NA_real_,
        signrep = NULL,
        trendfit = NULL,
        surrogate_method = surrogate_method,
        signal_mode = signal_mode
      ))
    }
    idx_active <- which(sig != 0)
  } else {
    sig_sd <- stats::sd(sig, na.rm = TRUE)
    if (!is.finite(sig_sd) || sig_sd < .Machine$double.eps) {
      if (warnings) warning("Signal has near-zero variance: surrogate distribution unavailable.")
      return(list(
        oscore_rp = NA_real_,
        fosc_rp = NA_real_,
        signrep = NULL,
        trendfit = NULL,
        surrogate_method = surrogate_method,
        signal_mode = signal_mode
      ))
    }
    idx_active <- which(abs(sig) > .Machine$double.eps)
  }
  if (length(idx_active) == 0) {
    if (warnings) warning("No active samples found for surrogate generation.")
    return(list(
      oscore_rp = NA_real_,
      fosc_rp = NA_real_,
      signrep = NULL,
      trendfit = NULL,
      surrogate_method = surrogate_method,
      signal_mode = signal_mode
    ))
  }

  tstart <- max(1L, idx_active[1] - 1L)
  tend <- min(length(sig), tail(idx_active, 1) + 1L)
  sigq <- sig[tstart:tend]

  trendfit <- list(distribution = "Shuffle", pd = NULL, pval = NA_real_, stats = NULL, trace = NULL)

  if (surrogate_method == "event_trend") {
    times <- which(sigq > 0)
    fits <- lapply(trend_dist, function(d) {
      tryCatch({
        fd <- fitdistrplus::fitdist(times, d)
        gs <- suppressWarnings(fitdistrplus::gofstat(fd))
        list(fd = fd, pval = gs$chisqpvalue)
      }, error = function(e) list(fd = NULL, pval = NA_real_))
    })
    pvals <- vapply(fits, function(f) f$pval, numeric(1))
    if (all(is.na(pvals))) {
      if (warnings) warning("Cannot fit distribution, shuffling data instead.")
      surrogate_method <- "event_shuffle"
    } else {
      id <- which.max(pvals)
      if (!is.na(pvals[id]) && pvals[id] >= trend_alpha) {
        fd <- fits[[id]]$fd
        xidx <- seq_len(length(sigq))
        pdf_vals <- fitdist_density(fd, xidx)
        pdf_vals[pdf_vals < 0] <- 0
        trace <- numeric(length(sig))
        trace[tstart:tend] <- pdf_vals / max(sum(pdf_vals), .Machine$double.eps)
        trendfit <- list(distribution = trend_dist[id], pd = fd, pval = pvals[id], stats = NULL, trace = trace)
      } else {
        if (warnings) warning("No suitable reference distribution found, using shuffle.")
        surrogate_method <- "event_shuffle"
      }
    }
  }

  oscore_rp <- rep(NA_real_, nrep)
  fosc_rp <- rep(NA_real_, nrep)
  signrep <- vector("list", nrep)

  for (rp in seq_len(nrep)) {
    if (surrogate_method == "phase_randomized") {
      sigrp <- phase_randomize_signal(sigq)
    } else if (surrogate_method == "event_trend" && !is.null(trendfit$pd)) {
      n_events <- max(0L, as.integer(round(sum(pmax(sigq, 0), na.rm = TRUE))))
      xidx <- seq_len(length(sigq))
      pdf_vals <- fitdist_density(trendfit$pd, xidx)
      pdf_vals[pdf_vals < 0] <- 0
      probs <- pdf_vals / max(sum(pdf_vals), .Machine$double.eps)
      sampled <- sample(xidx, size = n_events, replace = TRUE, prob = probs)
      sigrp <- tabulate(sampled, nbins = length(sigq))
    } else {
      sigrp <- numeric(length(sigq))
      events <- which(sigq > 0)
      if (length(events) > 0) {
        weights <- pmax(1L, as.integer(round(sigq[events])))
        ev_idx <- rep(events, weights)
        wind <- if (is.null(fpeak) || is.na(fpeak) || !is.finite(fpeak) || fpeak <= 0) fs / flim[1] else fs / fpeak
        for (i in ev_idx) {
          j <- round(stats::runif(1) * wind - wind / 2) + as.integer(i)
          j <- max(1L, min(length(sigrp), j))
          sigrp[j] <- sigrp[j] + 1
        }
      }
    }

    sig_full <- numeric(length(sig))
    sig_full[tstart:tend] <- sigrp
    res <- oscillation_score(sig_full,
                             fs = fs,
                             flim = flim,
                             fpeak = fpeak,
                             warnings = FALSE,
                             fcor = fcor,
                             taper = taper,
                             signal_mode = signal_mode)
    oscore_rp[rp] <- res$oscore
    fosc_rp[rp] <- res$fosc
    signrep[[rp]] <- sig_full
  }

  list(
    oscore_rp = oscore_rp,
    fosc_rp = fosc_rp,
    signrep = signrep,
    trendfit = trendfit,
    surrogate_method = surrogate_method,
    signal_mode = signal_mode
  )
}

phase_randomize_signal <- function(x) {
  sig <- as.numeric(x)
  sig[!is.finite(sig)] <- 0
  n <- length(sig)
  if (n <= 2) return(sig)

  mu <- mean(sig)
  centered <- sig - mu
  if (!any(abs(centered) > .Machine$double.eps)) {
    return(sig)
  }

  X <- fft(centered)
  if (n %% 2 == 0) {
    pos <- 2:(n / 2)
  } else {
    pos <- 2:((n + 1) / 2)
  }

  if (length(pos) > 0) {
    amp <- Mod(X[pos])
    phi <- stats::runif(length(pos), min = 0, max = 2 * pi)
    X_pos <- amp * exp(1i * phi)
    X[pos] <- X_pos
    X[n - pos + 2] <- Conj(X_pos)
  }

  X[1] <- Re(X[1])
  if (n %% 2 == 0) {
    X[n / 2 + 1] <- Re(X[n / 2 + 1])
  }

  as.numeric(Re(fft(X, inverse = TRUE) / n) + mu)
}

fitdist_density <- function(fd, x) {
  if (is.null(fd) || is.null(fd$estimate)) {
    stop("Invalid fitdist object: missing estimates.")
  }

  densfun <- NULL
  if (!is.null(fd$densfun) && (is.function(fd$densfun) || is.character(fd$densfun))) {
    densfun <- fd$densfun
  } else if (!is.null(fd$distname) && nzchar(fd$distname)) {
    cand <- paste0("d", fd$distname)
    if (exists(cand, mode = "function")) {
      densfun <- get(cand, mode = "function")
    }
  }
  if (is.null(densfun)) {
    stop("Could not resolve density function from fitdist object.")
  }

  do.call(densfun, c(list(x), as.list(fd$estimate)))
}
