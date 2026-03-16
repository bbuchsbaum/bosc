#' Grouped spectrum for binned subject-level series
#'
#' Computes per-subject spectra and their group average for a regular matrix of
#' subject-by-bin values or a grouped-series object.
#'
#' @param x numeric matrix with shape \code{subjects x bins}, or a list with
#'   fields \code{data}, \code{subjects}, and \code{bins}.
#' @param fs sampling rate in Hz along the bin axis.
#' @param flim optional length-2 frequency bounds.
#' @param detrend logical; if \code{TRUE}, remove a linear trend from each
#'   subject series before spectral estimation.
#' @param detrend_order non-negative integer polynomial order used when
#'   \code{detrend = TRUE}. The default \code{1} removes a linear trend, while
#'   higher orders remove higher-order polynomial structure.
#' @param pad_to optional FFT length used for zero-padding each subject series
#'   before spectral estimation.
#' @param spectrum spectrum scaling passed to \code{\link{spectral_peak}}.
#' @param taper taper passed to \code{\link{spectral_peak}}.
#' @param average logical; if \code{TRUE}, include the group-average spectrum.
#' @param fcor logical; apply 1/f correction to each subject spectrum.
#' @return list with \code{observed} and \code{meta} fields.
#' @export
group_spectrum <- function(x,
                           fs,
                           flim = NULL,
                           detrend = TRUE,
                           detrend_order = 1L,
                           pad_to = NULL,
                           spectrum = c("amplitude", "raw_power"),
                           taper = c("none", "hann", "hanning"),
                           average = TRUE,
                           fcor = FALSE) {
  spectrum <- match.arg(spectrum)
  taper <- match.arg(taper)
  gx <- as_grouped_series(x)
  validate_group_sampling(fs, flim, detrend_order)

  spec <- compute_group_spectrum(
    gx$data,
    fs = fs,
    flim = flim,
    detrend = detrend,
    detrend_order = detrend_order,
    pad_to = pad_to,
    spectrum = spectrum,
    taper = taper,
    fcor = fcor
  )

  observed <- list(
    freq = spec$freq,
    power = if (isTRUE(average)) spec$power else NULL,
    per_subject = spec$per_subject,
    peak = spec$peak,
    subject_peak = spec$subject_peak
  )

  list(
    observed = observed,
    meta = c(
      gx$meta,
      list(
        fs = fs,
        flim = flim,
        detrend = detrend,
        detrend_order = detrend_order,
        pad_to = pad_to,
        spectrum = spectrum,
        taper = taper,
      fcor = fcor,
      average = average,
      n_subjects = nrow(gx$data),
      n_bins = ncol(gx$data),
      measure = gx$measure,
      bins = gx$bins,
      subjects = gx$subjects
    )
  )
  )
}

#' Grouped spectrum with resampled null testing
#'
#' Builds a null distribution for the group-average spectrum by resampling each
#' subject series, then computes empirical p-values, z-scores, and optional FDR
#' correction across frequencies.
#'
#' @param x numeric matrix with shape \code{subjects x bins}, or a grouped-series
#'   object.
#' @param fs sampling rate in Hz along the bin axis.
#' @param flim optional length-2 frequency bounds.
#' @param null null generator; one of \code{"shuffle_labels"} or
#'   \code{"circular_shift"}.
#' @param nrep number of null draws.
#' @param detrend logical; if \code{TRUE}, remove a linear trend from each
#'   subject series before spectral estimation.
#' @param detrend_order non-negative integer polynomial order used when
#'   \code{detrend = TRUE}.
#' @param pad_to optional FFT length used for zero-padding each subject series
#'   before spectral estimation.
#' @param spectrum spectrum scaling passed to \code{\link{spectral_peak}}.
#' @param taper taper passed to \code{\link{spectral_peak}}.
#' @param p_adjust p-value adjustment; one of \code{"none"} or \code{"fdr"}.
#' @param fcor logical; apply 1/f correction to each subject spectrum.
#' @param seed optional random seed for reproducible null draws.
#' @return list with \code{observed}, \code{null}, \code{stats}, and \code{meta}.
#' @export
group_spectrum_test <- function(x,
                                fs,
                                flim = NULL,
                                null = c("shuffle_labels", "circular_shift"),
                                nrep = 1000,
                                detrend = TRUE,
                                detrend_order = 1L,
                                pad_to = NULL,
                                spectrum = c("amplitude", "raw_power"),
                                taper = c("none", "hann", "hanning"),
                                p_adjust = c("none", "fdr"),
                                fcor = FALSE,
                                seed = NULL) {
  null <- match.arg(null)
  spectrum <- match.arg(spectrum)
  taper <- match.arg(taper)
  p_adjust <- match.arg(p_adjust)
  if (length(nrep) != 1 || !is.finite(nrep) || nrep <= 0) {
    stop("nrep must be a positive scalar.")
  }
  nrep <- as.integer(round(nrep))
  if (!is.null(seed) && (length(seed) != 1 || !is.finite(seed))) {
    stop("seed must be NULL or a finite scalar.")
  }

  gx <- as_grouped_series(x)
  validate_group_sampling(fs, flim, detrend_order)
  obs <- compute_group_spectrum(
    gx$data,
    fs = fs,
    flim = flim,
    detrend = detrend,
    detrend_order = detrend_order,
    pad_to = pad_to,
    spectrum = spectrum,
    taper = taper,
    fcor = fcor
  )

  null_power <- with_seed(seed, {
    replicate(nrep, {
      draw <- generate_group_spectrum_null(gx$data, method = null)
      compute_group_spectrum(
        draw,
        fs = fs,
        flim = flim,
        detrend = detrend,
        detrend_order = detrend_order,
        pad_to = pad_to,
        spectrum = spectrum,
        taper = taper,
        fcor = fcor
      )$power
    })
  })
  if (is.null(dim(null_power))) {
    null_power <- matrix(null_power, nrow = length(obs$power), ncol = 1L)
  }

  p <- rowMeans(null_power >= obs$power)
  z <- (obs$power - rowMeans(null_power)) / apply(null_power, 1, stats::sd)
  z[!is.finite(z)] <- NA_real_
  p_adj <- if (p_adjust == "fdr") stats::p.adjust(p, method = "BH") else p

  list(
    observed = list(
      freq = obs$freq,
      power = obs$power,
      per_subject = obs$per_subject,
      peak = obs$peak,
      subject_peak = obs$subject_peak
    ),
    null = list(
      power = t(null_power),
      method = null,
      nrep = nrep
    ),
    stats = list(
      p = p,
      z = z,
      p_adj = p_adj,
      significant = p_adj <= 0.05
    ),
    meta = c(
      gx$meta,
      list(
        fs = fs,
        flim = flim,
        detrend = detrend,
        detrend_order = detrend_order,
        pad_to = pad_to,
        spectrum = spectrum,
        taper = taper,
        fcor = fcor,
        null = null,
        nrep = nrep,
        p_adjust = p_adjust,
        n_subjects = nrow(gx$data),
        n_bins = ncol(gx$data),
        measure = gx$measure,
        bins = gx$bins,
        subjects = gx$subjects
      )
    )
  )
}

#' Grouped spectrum with trial-level permutation nulls
#'
#' Starts from a long trial table, aggregates within bins to form per-subject
#' series, and builds a null distribution by permuting bin labels within each
#' subject or subject-by-condition group before re-aggregating.
#'
#' This workflow is intended for dense-sampling designs where trial timing is
#' varied across a regular grid and oscillatory structure is tested on the
#' resulting subject-by-bin series. The package implementation is generic, but
#' the memory-encoding application in Biba et al. (2026) is a concrete example
#' of this analysis pattern.
#'
#' @param data data.frame containing trial-wise values.
#' @param value character scalar naming the trial-level value column.
#' @param bin character scalar naming the bin column.
#' @param subject character scalar naming the subject column.
#' @param by optional character vector of additional grouping columns.
#' @param bins optional vector of target bins.
#' @param fun summary function used to aggregate trials within bins.
#' @param na_rm logical; if \code{TRUE}, remove missing trial values before
#'   aggregation.
#' @param complete logical; if \code{TRUE}, request explicit rows for missing
#'   group/bin combinations.
#' @param incomplete handling for incomplete observed grouped series; one of
#'   \code{"error"} or \code{"drop"}.
#' @param fs sampling rate in Hz along the bin axis.
#' @param flim optional length-2 frequency bounds.
#' @param null trial-level null generator; one of \code{"shuffle_bins"} or
#'   \code{"circular_shift_bins"}.
#' @param order_by optional within-group ordering column required for
#'   \code{"circular_shift_bins"}.
#' @param nrep number of null draws.
#' @param repair_incomplete how incomplete permuted groups are handled.
#'   \code{"none"} applies the same incomplete policy as the observed series.
#'   \code{"resample"} redraws each group's permutation until it yields a
#'   complete aggregated series.
#' @param max_resample maximum redraw attempts per group when
#'   \code{repair_incomplete = "resample"}.
#' @param detrend logical; if \code{TRUE}, detrend each series before spectral
#'   estimation.
#' @param detrend_order polynomial detrend order used when \code{detrend = TRUE}.
#' @param pad_to optional FFT length used for zero-padding each subject series
#'   before spectral estimation.
#' @param spectrum spectrum scaling passed to \code{\link{spectral_peak}}.
#' @param taper taper passed to \code{\link{spectral_peak}}.
#' @param p_adjust p-value adjustment; one of \code{"none"} or \code{"fdr"}.
#' @param fcor logical; apply 1/f correction to each subject spectrum.
#' @param seed optional random seed.
#' @return list with \code{observed}, \code{null}, \code{stats}, and
#'   \code{meta}. The \code{meta} field includes the aggregated observed
#'   grouped-series object.
#' @references
#' Biba, T. M., Decker, A., Herrmann, B., Fukuda, K., Katz, C. N.,
#' Valiante, T. A., & Duncan, K. (2026). Episodic memory encoding
#' fluctuates at a theta rhythm of 3-10 Hz. *Nature Human Behaviour*.
#' \doi{10.1038/s41562-026-02416-5}
#' @keywords internal
group_spectrum_test_trials <- function(data,
                                       value,
                                       bin,
                                       subject,
                                       by = NULL,
                                       bins = NULL,
                                       fun = mean,
                                       na_rm = TRUE,
                                       complete = TRUE,
                                       incomplete = c("error", "drop"),
                                       fs,
                                       flim = NULL,
                                       null = c("shuffle_bins", "circular_shift_bins"),
                                       order_by = NULL,
                                       nrep = 1000,
                                       repair_incomplete = c("none", "resample"),
                                       max_resample = 100L,
                                       detrend = TRUE,
                                       detrend_order = 1L,
                                       pad_to = NULL,
                                       spectrum = c("amplitude", "raw_power"),
                                       taper = c("none", "hann", "hanning"),
                                       p_adjust = c("none", "fdr"),
                                       fcor = FALSE,
                                       seed = NULL) {
  data <- validate_data_frame(data)
  value <- validate_column_name(value, "value")
  bin <- validate_column_name(bin, "bin")
  subject <- validate_column_name(subject, "subject")
  by <- validate_optional_columns(by, data, "by")
  if (!value %in% names(data)) stop("value column not found in data.")
  if (!bin %in% names(data)) stop("bin column not found in data.")
  if (!subject %in% names(data)) stop("subject column not found in data.")
  null <- match.arg(null)
  repair_incomplete <- match.arg(repair_incomplete)
  spectrum <- match.arg(spectrum)
  taper <- match.arg(taper)
  p_adjust <- match.arg(p_adjust)
  incomplete <- match.arg(incomplete)
  if (length(nrep) != 1 || !is.finite(nrep) || nrep <= 0) {
    stop("nrep must be a positive scalar.")
  }
  nrep <- as.integer(round(nrep))
  if (length(max_resample) != 1 || !is.finite(max_resample) || max_resample <= 0) {
    stop("max_resample must be a positive scalar.")
  }
  max_resample <- as.integer(round(max_resample))

  validate_group_sampling(fs, flim, detrend_order)
  group_cols <- c(subject, by)
  observed_series <- make_grouped_series(
    data = data,
    value = value,
    bin = bin,
    subject = subject,
    by = by,
    bins = bins,
    fun = fun,
    na_rm = na_rm,
    complete = complete,
    incomplete = incomplete
  )

  trial_data <- data
  if (incomplete == "drop") {
    trial_data <- filter_to_group_keys(trial_data, observed_series$group_keys, group_cols)
  }

  obs <- compute_group_spectrum(
    observed_series$data,
    fs = fs,
    flim = flim,
    detrend = detrend,
    detrend_order = detrend_order,
    pad_to = pad_to,
    spectrum = spectrum,
    taper = taper,
    fcor = fcor
  )

  perm_method <- if (null == "shuffle_bins") "shuffle" else "circular_shift"
  if (perm_method == "circular_shift" && is.null(order_by)) {
    stop("order_by must be provided for circular_shift_bins.")
  }

  null_power <- with_seed(seed, {
    replicate(nrep, {
      if (repair_incomplete == "resample") {
        permuted <- permute_complete_trial_groups(
          data = trial_data,
          column = bin,
          group_cols = group_cols,
          group_keys = observed_series$group_keys,
          value = value,
          bins = observed_series$bins,
          fun = fun,
          na_rm = na_rm,
          out = observed_series$measure,
          method = perm_method,
          order_by = order_by,
          max_resample = max_resample
        )
        series <- make_grouped_series(
          data = permuted,
          value = value,
          bin = bin,
          subject = subject,
          by = by,
          bins = observed_series$bins,
          fun = fun,
          na_rm = na_rm,
          complete = complete,
          incomplete = "error",
          out = observed_series$measure
        )
      } else {
        permuted <- permute_trial_labels(
          data = trial_data,
          column = bin,
          by = group_cols,
          method = perm_method,
          order_by = order_by
        )
        series <- make_grouped_series(
          data = permuted,
          value = value,
          bin = bin,
          subject = subject,
          by = by,
          bins = observed_series$bins,
          fun = fun,
          na_rm = na_rm,
          complete = complete,
          incomplete = incomplete,
          out = observed_series$measure
        )
      }
      compute_group_spectrum(
        series$data,
        fs = fs,
        flim = flim,
        detrend = detrend,
        detrend_order = detrend_order,
        pad_to = pad_to,
        spectrum = spectrum,
        taper = taper,
        fcor = fcor
      )$power
    })
  })

  if (is.null(dim(null_power))) {
    null_power <- matrix(null_power, nrow = length(obs$power), ncol = 1L)
  }

  p <- rowMeans(null_power >= obs$power)
  z <- (obs$power - rowMeans(null_power)) / apply(null_power, 1, stats::sd)
  z[!is.finite(z)] <- NA_real_
  p_adj <- if (p_adjust == "fdr") stats::p.adjust(p, method = "BH") else p

  list(
    observed = list(
      freq = obs$freq,
      power = obs$power,
      per_subject = obs$per_subject,
      peak = obs$peak,
      subject_peak = obs$subject_peak
    ),
    null = list(
      power = t(null_power),
      method = null,
      nrep = nrep
    ),
    stats = list(
      p = p,
      z = z,
      p_adj = p_adj,
      significant = p_adj <= 0.05
    ),
    meta = list(
      fs = fs,
      flim = flim,
      detrend = detrend,
      detrend_order = detrend_order,
      pad_to = pad_to,
      spectrum = spectrum,
      taper = taper,
      fcor = fcor,
      null = null,
      nrep = nrep,
      p_adjust = p_adjust,
      value = value,
      bin = bin,
      subject = subject,
      by = by,
      incomplete = incomplete,
      repair_incomplete = repair_incomplete,
      grouped_series = observed_series
    )
  )
}

as_grouped_series <- function(x) {
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("x matrix must be numeric.")
    if (nrow(x) < 1 || ncol(x) < 2) {
      stop("x must have at least one subject and two bins.")
    }
    if (any(!is.finite(x))) {
      stop("x must contain only finite values.")
    }
    return(list(
      data = unname(x),
      subjects = if (is.null(rownames(x))) as.character(seq_len(nrow(x))) else rownames(x),
      bins = seq_len(ncol(x)),
      measure = "value",
      meta = list(input = "matrix")
    ))
  }

  if (is.list(x) && !is.null(x$data)) {
    dat <- x$data
    if (!is.matrix(dat) || !is.numeric(dat)) {
      stop("x$data must be a numeric matrix.")
    }
    if (nrow(dat) < 1 || ncol(dat) < 2) {
      stop("x$data must have at least one subject and two bins.")
    }
    if (any(!is.finite(dat))) {
      stop("x$data must contain only finite values.")
    }
    subjects <- x$subjects
    if (is.null(subjects)) {
      subjects <- if (is.null(rownames(dat))) as.character(seq_len(nrow(dat))) else rownames(dat)
    }
    bins <- x$bins
    if (is.null(bins)) bins <- seq_len(ncol(dat))
    if (length(subjects) != nrow(dat)) {
      stop("Length of x$subjects must match nrow(x$data).")
    }
    if (length(bins) != ncol(dat)) {
      stop("Length of x$bins must match ncol(x$data).")
    }
    return(list(
      data = unname(dat),
      subjects = as.character(subjects),
      bins = bins,
      measure = if (!is.null(x$measure)) x$measure else "value",
      meta = if (!is.null(x$meta)) x$meta else list(input = "grouped_series")
    ))
  }

  stop("x must be a numeric matrix or grouped-series list.")
}

validate_group_sampling <- function(fs, flim, detrend_order = 1L) {
  if (length(fs) != 1 || !is.finite(fs) || fs <= 0) {
    stop("fs must be a positive scalar.")
  }
  if (!is.null(flim)) {
    if (!is.numeric(flim) || length(flim) != 2 || any(!is.finite(flim))) {
      stop("flim must be NULL or a length-2 finite numeric vector.")
    }
    if (flim[1] <= 0 || flim[2] <= 0 || flim[1] >= flim[2]) {
      stop("flim must contain increasing positive frequencies.")
    }
  }
  if (length(detrend_order) != 1 || !is.finite(detrend_order) || detrend_order < 0) {
    stop("detrend_order must be a non-negative scalar.")
  }
}

compute_group_spectrum <- function(x, fs, flim, detrend, detrend_order, pad_to = NULL, spectrum = "amplitude", taper, fcor) {
  per_subject <- vector("list", nrow(x))
  freq <- NULL
  for (i in seq_len(nrow(x))) {
    xi <- preprocess_group_series(x[i, ], detrend = detrend, detrend_order = detrend_order)
    sp <- spectral_peak(
      xi,
      fs = fs,
      flim = flim,
      fcor = fcor,
      pad_to = pad_to,
      spectrum = spectrum,
      taper = taper
    )
    if (length(sp$fxx) == 0 || length(sp$spectrum) == 0) {
      stop("No spectral bins available for the requested settings.")
    }
    if (is.null(freq)) {
      freq <- sp$fxx
    } else if (!isTRUE(all.equal(freq, sp$fxx, tolerance = 1e-12))) {
      stop("Incompatible frequency grids across subjects.")
    }
    per_subject[[i]] <- sp
  }

  power_mat <- do.call(rbind, lapply(per_subject, function(sp) sp$spectrum))
  subject_peak <- vapply(per_subject, function(sp) sp$freq, numeric(1))
  mean_power <- colMeans(power_mat)
  peak <- freq[which.max(mean_power)]

  list(
    freq = freq,
    power = mean_power,
    per_subject = power_mat,
    peak = peak,
    subject_peak = subject_peak
  )
}

preprocess_group_series <- function(x, detrend = TRUE, detrend_order = 1L) {
  x <- as.numeric(x)
  if (!isTRUE(detrend)) return(x)
  t <- seq_along(x)
  detrend_order <- as.integer(round(detrend_order))
  form <- if (detrend_order == 0L) {
    stats::as.formula("x ~ 1")
  } else {
    stats::as.formula(
      paste0("x ~ poly(t, ", detrend_order, ", raw = TRUE)")
    )
  }
  fit <- stats::lm(form, data = data.frame(x = x, t = t))
  stats::residuals(fit)
}

generate_group_spectrum_null <- function(x, method = c("shuffle_labels", "circular_shift")) {
  method <- match.arg(method)
  out <- x
  for (i in seq_len(nrow(x))) {
    if (method == "shuffle_labels") {
      out[i, ] <- sample(x[i, ], length(x[i, ]), replace = FALSE)
    } else if (method == "circular_shift") {
      k <- sample.int(ncol(x), size = 1) - 1L
      out[i, ] <- circ_shift_vector(x[i, ], k)
    }
  }
  out
}

circ_shift_vector <- function(x, k) {
  n <- length(x)
  if (n == 0) return(x)
  k <- k %% n
  if (k == 0) return(x)
  c(x[(n - k + 1):n], x[seq_len(n - k)])
}

filter_to_group_keys <- function(data, keys, cols) {
  keep <- rep(FALSE, nrow(data))
  for (i in seq_len(nrow(keys))) {
    mask <- rep(TRUE, nrow(data))
    for (col in cols) {
      mask <- mask & same_or_na(data[[col]], keys[[col]][i])
    }
    keep <- keep | mask
  }
  data[keep, , drop = FALSE]
}

permute_complete_trial_groups <- function(data,
                                          column,
                                          group_cols,
                                          group_keys,
                                          value,
                                          bins,
                                          fun,
                                          na_rm,
                                          out,
                                          method,
                                          order_by,
                                          max_resample) {
  out_groups <- vector("list", nrow(group_keys))
  for (i in seq_len(nrow(group_keys))) {
    group_data <- subset_group_rows(data, group_keys, i, group_cols)
    success <- FALSE
    for (attempt in seq_len(max_resample)) {
      permuted <- permute_trial_labels(
        data = group_data,
        column = column,
        by = NULL,
        method = method,
        order_by = order_by
      )
      agg <- aggregate_by_bin(
        data = permuted,
        value = value,
        bin = column,
        by = NULL,
        bins = bins,
        fun = fun,
        na_rm = na_rm,
        complete = TRUE,
        out = out
      )
      if (nrow(agg) == length(bins) && all(is.finite(agg[[out]]))) {
        out_groups[[i]] <- permuted
        success <- TRUE
        break
      }
    }
    if (!success) {
      stop("Could not generate a complete permuted series within max_resample.")
    }
  }
  do.call(rbind, out_groups)
}

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(eval.parent(substitute(expr)))
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  eval.parent(substitute(expr))
}
