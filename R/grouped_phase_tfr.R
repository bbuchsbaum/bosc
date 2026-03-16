#' Grouped time-frequency map via narrowband Hilbert decomposition
#'
#' Computes a grouped time-frequency amplitude map by filtering each subject
#' series at a sequence of center frequencies and averaging the resulting
#' narrowband amplitudes.
#'
#' @param x numeric matrix with shape \code{subjects x bins}, or a grouped-series
#'   object.
#' @param fs sampling rate in Hz along the bin axis.
#' @param freqs numeric vector of center frequencies.
#' @param bandwidth optional half-bandwidth in Hz. If \code{NULL}, a default is
#'   derived from \code{freqs}.
#' @param detrend logical; if \code{TRUE}, detrend each subject series first.
#' @param reflect logical; if \code{TRUE}, reflect the series at both ends before
#'   filtering and trim back to the original span.
#' @param edge number of bins to trim from each edge after filtering.
#' @param filtorder,demean,tol passed to \code{\link{narrowband_hilbert}}.
#' @return list with \code{observed} and \code{meta}.
#' @export
group_tfr <- function(x,
                      fs,
                      freqs,
                      bandwidth = NULL,
                      detrend = TRUE,
                      reflect = FALSE,
                      edge = 0,
                      filtorder = 2,
                      demean = TRUE,
                      tol = 1e2) {
  gx <- as_grouped_series(x)
  freqs <- validate_group_freqs(freqs, fs)
  bandwidth <- resolve_group_bandwidth(freqs, bandwidth, fs)
  edge <- validate_group_edge(edge, ncol(gx$data))

  tf <- compute_group_tfr(
    gx$data,
    fs = fs,
    freqs = freqs,
    bandwidth = bandwidth,
    detrend = detrend,
    reflect = reflect,
    edge = edge,
    filtorder = filtorder,
    demean = demean,
    tol = tol
  )

  list(
    observed = list(
      map = tf$map,
      per_subject = tf$per_subject,
      time = gx$bins[tf$keep_idx],
      freq = freqs
    ),
    meta = c(
      gx$meta,
      list(
        fs = fs,
        freqs = freqs,
        bandwidth = bandwidth,
        detrend = detrend,
        reflect = reflect,
        edge = edge,
        filtorder = filtorder,
        demean = demean,
        tol = tol,
        measure = gx$measure,
        bins = gx$bins[tf$keep_idx],
        subjects = gx$subjects,
        n_subjects = nrow(gx$data),
        n_bins = length(tf$keep_idx)
      )
    )
  )
}

#' Grouped time-frequency map with resampled null testing
#'
#' Computes an empirical null distribution for grouped narrowband amplitude maps
#' and returns cell-wise statistics with optional cluster summaries.
#'
#' @param x numeric matrix with shape \code{subjects x bins}, or a grouped-series
#'   object.
#' @param fs sampling rate in Hz along the bin axis.
#' @param freqs numeric vector of center frequencies.
#' @param bandwidth optional half-bandwidth in Hz.
#' @param null null generator; one of \code{"shuffle_labels"} or
#'   \code{"circular_shift"}.
#' @param nrep number of null draws.
#' @param detrend,reflect,edge,filtorder,demean,tol passed to
#'   \code{\link{group_tfr}}.
#' @param p_adjust p-value adjustment; one of \code{"none"} or \code{"fdr"}.
#' @param cluster logical; if \code{TRUE}, detect suprathreshold clusters in the
#'   z-map.
#' @param cluster_threshold z-threshold used for cluster detection.
#' @param seed optional random seed.
#' @return list with \code{observed}, \code{null}, \code{stats}, and \code{meta}.
#' @export
group_tfr_test <- function(x,
                           fs,
                           freqs,
                           bandwidth = NULL,
                           null = c("shuffle_labels", "circular_shift"),
                           nrep = 1000,
                           detrend = TRUE,
                           reflect = FALSE,
                           edge = 0,
                           filtorder = 2,
                           demean = TRUE,
                           tol = 1e2,
                           p_adjust = c("none", "fdr"),
                           cluster = TRUE,
                           cluster_threshold = stats::qnorm(0.975),
                           seed = NULL) {
  null <- match.arg(null)
  p_adjust <- match.arg(p_adjust)
  if (length(nrep) != 1 || !is.finite(nrep) || nrep <= 0) {
    stop("nrep must be a positive scalar.")
  }
  if (length(cluster_threshold) != 1 || !is.finite(cluster_threshold) || cluster_threshold <= 0) {
    stop("cluster_threshold must be a positive scalar.")
  }
  nrep <- as.integer(round(nrep))

  gx <- as_grouped_series(x)
  freqs <- validate_group_freqs(freqs, fs)
  bandwidth <- resolve_group_bandwidth(freqs, bandwidth, fs)
  edge <- validate_group_edge(edge, ncol(gx$data))

  obs <- compute_group_tfr(
    gx$data,
    fs = fs,
    freqs = freqs,
    bandwidth = bandwidth,
    detrend = detrend,
    reflect = reflect,
    edge = edge,
    filtorder = filtorder,
    demean = demean,
    tol = tol
  )

  null_maps <- with_seed(seed, {
    arr <- array(
      NA_real_,
      dim = c(nrep, length(freqs), length(obs$keep_idx))
    )
    for (i in seq_len(nrep)) {
      draw <- generate_group_spectrum_null(gx$data, method = null)
      arr[i, , ] <- compute_group_tfr(
        draw,
        fs = fs,
        freqs = freqs,
        bandwidth = bandwidth,
        detrend = detrend,
        reflect = reflect,
        edge = edge,
        filtorder = filtorder,
        demean = demean,
        tol = tol
      )$map
    }
    arr
  })

  null_mean <- apply(null_maps, c(2, 3), mean)
  null_sd <- apply(null_maps, c(2, 3), stats::sd)
  z <- (obs$map - null_mean) / null_sd
  z[!is.finite(z)] <- NA_real_

  p <- matrix(NA_real_, nrow = nrow(obs$map), ncol = ncol(obs$map))
  for (fi in seq_len(nrow(obs$map))) {
    for (ti in seq_len(ncol(obs$map))) {
      p[fi, ti] <- mean(null_maps[, fi, ti] >= obs$map[fi, ti])
    }
  }
  p_adj <- if (p_adjust == "fdr") {
    matrix(stats::p.adjust(as.vector(p), method = "BH"), nrow = nrow(p), ncol = ncol(p))
  } else {
    p
  }

  clusters <- list()
  if (isTRUE(cluster)) {
    clusters <- detect_clusters(
      z,
      threshold = cluster_threshold,
      method = "extract",
      freq_axis = freqs,
      time_axis = gx$bins[obs$keep_idx]
    )
  }

  list(
    observed = list(
      map = obs$map,
      per_subject = obs$per_subject,
      time = gx$bins[obs$keep_idx],
      freq = freqs
    ),
    null = list(
      map = null_maps,
      method = null,
      nrep = nrep
    ),
    stats = list(
      z = z,
      p = p,
      p_adj = p_adj,
      significant = p_adj <= 0.05,
      clusters = clusters
    ),
    meta = c(
      gx$meta,
      list(
        fs = fs,
        freqs = freqs,
        bandwidth = bandwidth,
        detrend = detrend,
        reflect = reflect,
        edge = edge,
        filtorder = filtorder,
        demean = demean,
        tol = tol,
        null = null,
        nrep = nrep,
        p_adjust = p_adjust,
        cluster = cluster,
        cluster_threshold = cluster_threshold,
        measure = gx$measure,
        bins = gx$bins[obs$keep_idx],
        subjects = gx$subjects,
        n_subjects = nrow(gx$data),
        n_bins = length(obs$keep_idx)
      )
    )
  )
}

#' Grouped phase consistency over bins
#'
#' Computes per-bin phase concentration across subjects after narrowband Hilbert
#' decomposition of subject-level series.
#'
#' @param x numeric matrix with shape \code{subjects x bins}, or a grouped-series
#'   object.
#' @param fs sampling rate in Hz along the bin axis.
#' @param freq center frequency in Hz.
#' @param bandwidth half-bandwidth in Hz.
#' @param detrend,reflect,edge,filtorder,demean,tol preprocessing controls.
#' @return list with \code{observed}, \code{stats}, and \code{meta}.
#' @export
group_phase_consistency <- function(x,
                                    fs,
                                    freq,
                                    bandwidth = 0.5,
                                    detrend = TRUE,
                                    reflect = FALSE,
                                    edge = 0,
                                    filtorder = 2,
                                    demean = TRUE,
                                    tol = 1e2) {
  gx <- as_grouped_series(x)
  validate_group_freqs(freq, fs)
  edge <- validate_group_edge(edge, ncol(gx$data))

  ph <- compute_group_phase_matrix(
    gx$data,
    fs = fs,
    freq = freq,
    bandwidth = bandwidth,
    detrend = detrend,
    reflect = reflect,
    edge = edge,
    filtorder = filtorder,
    demean = demean,
    tol = tol
  )

  r <- apply(ph$phase, 2, circ_r)
  mu <- apply(ph$phase, 2, circ_mean)
  p <- vapply(seq_len(ncol(ph$phase)), function(i) {
    circ_rayleigh(ph$phase[, i])$pval
  }, numeric(1))

  list(
    observed = list(
      phase = ph$phase,
      mean_phase = mu,
      time = gx$bins[ph$keep_idx]
    ),
    stats = list(
      r = r,
      p = p,
      significant = p <= 0.05
    ),
    meta = c(
      gx$meta,
      list(
        fs = fs,
        freq = freq,
        bandwidth = bandwidth,
        detrend = detrend,
        reflect = reflect,
        edge = edge,
        filtorder = filtorder,
        demean = demean,
        tol = tol,
        measure = gx$measure,
        bins = gx$bins[ph$keep_idx],
        subjects = gx$subjects,
        n_subjects = nrow(gx$data),
        n_bins = length(ph$keep_idx)
      )
    )
  )
}

#' Grouped phase difference over bins
#'
#' Computes per-bin phase differences between two matched grouped series objects
#' or matrices, then summarizes phase concentration of those differences across
#' subjects.
#'
#' @param x1,x2 numeric matrices with identical shape, or grouped-series objects
#'   with matching \code{subjects} and \code{bins}.
#' @param fs sampling rate in Hz along the bin axis.
#' @param freq center frequency in Hz.
#' @param bandwidth half-bandwidth in Hz.
#' @param mode one of \code{"signed_pi"} or \code{"wrapped_2pi"}.
#' @param detrend,reflect,edge,filtorder,demean,tol preprocessing controls.
#' @return list with \code{observed}, \code{stats}, and \code{meta}.
#' @export
group_phase_difference <- function(x1,
                                   x2,
                                   fs,
                                   freq,
                                   bandwidth = 0.5,
                                   mode = c("signed_pi", "wrapped_2pi"),
                                   detrend = TRUE,
                                   reflect = FALSE,
                                   edge = 0,
                                   filtorder = 2,
                                   demean = TRUE,
                                   tol = 1e2) {
  mode <- match.arg(mode)
  g1 <- as_grouped_series(x1)
  g2 <- as_grouped_series(x2)
  if (!identical(g1$subjects, g2$subjects)) {
    stop("x1 and x2 must have matching subject order.")
  }
  if (!isTRUE(all.equal(g1$bins, g2$bins, tolerance = 1e-12))) {
    stop("x1 and x2 must have matching bins.")
  }
  validate_group_freqs(freq, fs)
  edge <- validate_group_edge(edge, ncol(g1$data))

  ph1 <- compute_group_phase_matrix(
    g1$data,
    fs = fs,
    freq = freq,
    bandwidth = bandwidth,
    detrend = detrend,
    reflect = reflect,
    edge = edge,
    filtorder = filtorder,
    demean = demean,
    tol = tol
  )
  ph2 <- compute_group_phase_matrix(
    g2$data,
    fs = fs,
    freq = freq,
    bandwidth = bandwidth,
    detrend = detrend,
    reflect = reflect,
    edge = edge,
    filtorder = filtorder,
    demean = demean,
    tol = tol
  )

  diff <- wrap_phase(ph1$phase - ph2$phase, mode = mode)
  r <- apply(diff, 2, circ_r)
  mu <- apply(diff, 2, circ_mean)
  p <- vapply(seq_len(ncol(diff)), function(i) {
    circ_rayleigh(diff[, i])$pval
  }, numeric(1))

  list(
    observed = list(
      phase_difference = diff,
      mean_phase_difference = mu,
      time = g1$bins[ph1$keep_idx]
    ),
    stats = list(
      r = r,
      p = p,
      significant = p <= 0.05
    ),
    meta = c(
      g1$meta,
      list(
        fs = fs,
        freq = freq,
        bandwidth = bandwidth,
        mode = mode,
        detrend = detrend,
        reflect = reflect,
        edge = edge,
        filtorder = filtorder,
        demean = demean,
        tol = tol,
        bins = g1$bins[ph1$keep_idx],
        subjects = g1$subjects,
        n_subjects = nrow(g1$data),
        n_bins = length(ph1$keep_idx)
      )
    )
  )
}

compute_group_tfr <- function(x,
                              fs,
                              freqs,
                              bandwidth,
                              detrend,
                              reflect,
                              edge,
                              filtorder,
                              demean,
                              tol) {
  keep_idx <- seq_len(ncol(x))
  if (edge > 0) {
    keep_idx <- keep_idx[(edge + 1):(length(keep_idx) - edge)]
  }

  arr <- array(NA_real_, dim = c(nrow(x), length(freqs), length(keep_idx)))
  for (si in seq_len(nrow(x))) {
    xi <- preprocess_group_series(x[si, ], detrend = detrend)
    xf <- maybe_reflect_series(xi, reflect = reflect)
    trim_idx <- reflected_trim_idx(length(xi), reflect = reflect)
    for (fi in seq_along(freqs)) {
      band <- frequency_band(freqs[fi], bandwidth[fi], fs)
      nb <- narrowband_hilbert(
        data = xf,
        fs = fs,
        freqlim = band,
        tol = tol,
        filtorder = filtorder,
        demean = demean
      )
      amp <- Mod(nb$analytic)[trim_idx]
      arr[si, fi, ] <- amp[keep_idx]
    }
  }

  list(
    map = apply(arr, c(2, 3), mean, na.rm = TRUE),
    per_subject = arr,
    keep_idx = keep_idx
  )
}

compute_group_phase_matrix <- function(x,
                                       fs,
                                       freq,
                                       bandwidth,
                                       detrend,
                                       reflect,
                                       edge,
                                       filtorder,
                                       demean,
                                       tol) {
  keep_idx <- seq_len(ncol(x))
  if (edge > 0) {
    keep_idx <- keep_idx[(edge + 1):(length(keep_idx) - edge)]
  }

  out <- matrix(NA_real_, nrow = nrow(x), ncol = length(keep_idx))
  for (si in seq_len(nrow(x))) {
    xi <- preprocess_group_series(x[si, ], detrend = detrend)
    xf <- maybe_reflect_series(xi, reflect = reflect)
    trim_idx <- reflected_trim_idx(length(xi), reflect = reflect)
    band <- frequency_band(freq, bandwidth, fs)
    nb <- narrowband_hilbert(
      data = xf,
      fs = fs,
      freqlim = band,
      tol = tol,
      filtorder = filtorder,
      demean = demean
    )
    ph <- Arg(nb$analytic)[trim_idx]
    out[si, ] <- ph[keep_idx]
  }
  list(phase = out, keep_idx = keep_idx)
}

validate_group_freqs <- function(freqs, fs) {
  if (!is.numeric(freqs) || length(freqs) == 0 || any(!is.finite(freqs))) {
    stop("freqs must be a non-empty finite numeric vector.")
  }
  if (any(freqs <= 0) || any(freqs >= fs / 2)) {
    stop("freqs must lie strictly between 0 and Nyquist.")
  }
  as.numeric(freqs)
}

resolve_group_bandwidth <- function(freqs, bandwidth, fs) {
  if (is.null(bandwidth)) {
    if (length(freqs) >= 2) {
      bw <- rep(max(stats::median(diff(sort(unique(freqs)))) / 2, 0.25), length(freqs))
    } else {
      bw <- rep(max(freqs / 8, 0.5), length(freqs))
    }
  } else if (length(bandwidth) == 1) {
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      stop("bandwidth must be positive.")
    }
    bw <- rep(as.numeric(bandwidth), length(freqs))
  } else if (length(bandwidth) == length(freqs)) {
    if (any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
      stop("bandwidth values must be positive.")
    }
    bw <- as.numeric(bandwidth)
  } else {
    stop("bandwidth must be NULL, length 1, or the same length as freqs.")
  }
  if (any(freqs - bw <= 0) || any(freqs + bw >= fs / 2)) {
    stop("bandwidth produces invalid bands relative to Nyquist or zero.")
  }
  bw
}

validate_group_edge <- function(edge, n_bin) {
  if (length(edge) != 1 || !is.finite(edge) || edge < 0 || edge != floor(edge)) {
    stop("edge must be a non-negative integer.")
  }
  if (2 * edge >= n_bin) {
    stop("edge trims away all bins.")
  }
  as.integer(edge)
}

maybe_reflect_series <- function(x, reflect = FALSE) {
  if (!isTRUE(reflect) || length(x) < 2) return(x)
  c(rev(x[-1]), x, rev(x[-length(x)]))
}

reflected_trim_idx <- function(n, reflect = FALSE) {
  if (!isTRUE(reflect) || n < 2) return(seq_len(n))
  n:(2 * n - 1)
}

frequency_band <- function(freq, bandwidth, fs) {
  c(max(.Machine$double.eps, freq - bandwidth), min(fs / 2 - .Machine$double.eps, freq + bandwidth))
}

wrap_phase <- function(x, mode = c("signed_pi", "wrapped_2pi")) {
  mode <- match.arg(mode)
  if (mode == "signed_pi") {
    (x + pi) %% (2 * pi) - pi
  } else {
    x %% (2 * pi)
  }
}
