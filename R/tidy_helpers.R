#' Tidy summary of oscillation score results
#'
#' Converts an \code{\link{oscillation_score}} or \code{\link{oscillation_score_z}}
#' result into a one-row data.frame. If a list of results is provided, returns one
#' row per element.
#'
#' @param x a result list (or list of result lists).
#' @return data.frame with summary columns.
#' @export
oscore_tidy <- function(x) {
  if (is.null(x)) stop("x cannot be NULL.")

  is_result <- function(y) is.list(y) && !is.null(y$oscore) && !is.null(y$fosc)

  if (is_result(x)) {
    flim <- y <- x
    fmin <- if (!is.null(flim$flim) && length(flim$flim) == 2) flim$flim[1] else NA_real_
    fmax <- if (!is.null(flim$flim) && length(flim$flim) == 2) flim$flim[2] else NA_real_
    out <- list(
      oscore = y$oscore,
      fosc = y$fosc,
      fmin = fmin,
      fmax = fmax
    )
    if (!is.null(y$z)) out$z <- y$z
    if (!is.null(y$pval)) out$pval <- y$pval
    if (!is.null(y$significant)) out$significant <- y$significant
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }

  if (is.list(x) && length(x) > 0 && all(vapply(x, is_result, logical(1)))) {
    dfs <- lapply(x, oscore_tidy)
    out <- do.call(rbind, dfs)
    rownames(out) <- NULL
    return(out)
  }

  stop("x must be an oscillation score result or list of results.")
}

#' Tidy spectrum from oscillation score or spectral peak output
#'
#' Converts a result containing a spectrum into a long data.frame with columns
#' \code{freq} and \code{power}. Supports outputs from
#' \code{\link{oscillation_score}} (fields \code{freqs}/\code{flimfft}) and
#' \code{\link{spectral_peak}} (fields \code{fxx}/\code{spectrum}).
#'
#' If a list of results is provided, an \code{id} column is added.
#'
#' @param x a spectrum-containing result list (or list of such results).
#' @return data.frame with columns \code{freq}, \code{power}, and optionally \code{id}.
#' @export
oscore_spectrum <- function(x) {
  if (is.null(x)) stop("x cannot be NULL.")

  extract_one <- function(y) {
    if (is.null(y)) return(NULL)
    if (!is.null(y$freqs) && !is.null(y$flimfft)) {
      return(data.frame(freq = y$freqs, power = y$flimfft, stringsAsFactors = FALSE))
    }
    if (!is.null(y$fxx) && !is.null(y$spectrum)) {
      return(data.frame(freq = y$fxx, power = y$spectrum, stringsAsFactors = FALSE))
    }
    NULL
  }

  one <- extract_one(x)
  if (!is.null(one)) return(one)

  if (is.list(x) && length(x) > 0) {
    dfs <- lapply(seq_along(x), function(i) {
      df <- extract_one(x[[i]])
      if (is.null(df)) return(NULL)
      df$id <- i
      df
    })
    dfs <- Filter(Negate(is.null), dfs)
    if (length(dfs) == 0) stop("No spectrum found in x.")
    out <- do.call(rbind, dfs)
    rownames(out) <- NULL
    return(out)
  }

  stop("x must contain spectrum fields.")
}

#' Tidy summary of detected clusters
#'
#' Converts the list returned by \code{\link{detect_clusters}} into a data.frame
#' with one row per cluster.
#'
#' @param clusters list of clusters as returned by \code{\link{detect_clusters}}.
#' @return data.frame with cluster summaries.
#' @export
clusters_tidy <- function(clusters) {
  if (is.null(clusters) || length(clusters) == 0) {
    return(data.frame(
      cluster = integer(0),
      blobID = numeric(0),
      sign = integer(0),
      n_pixels = integer(0),
      time_min = numeric(0),
      time_max = numeric(0),
      freq_min = numeric(0),
      freq_max = numeric(0),
      CoMtime = numeric(0),
      CoMfreq = numeric(0),
      peakZscore = numeric(0),
      sumZscore = numeric(0),
      peakpAdj = numeric(0),
      minpAdj = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(seq_along(clusters), function(i) {
    cl <- clusters[[i]]
    pAdj <- cl$pAdj
    if (is.null(pAdj)) pAdj <- rep(NA_real_, length(cl$Zscores))
    pAdj_valid <- pAdj[is.finite(pAdj)]
    data.frame(
      cluster = i,
      blobID = cl$blobID,
      sign = if (!is.null(cl$blobID)) sign(cl$blobID) else sign(cl$sumZscore),
      n_pixels = length(cl$Zscores),
      time_min = min(cl$times, na.rm = TRUE),
      time_max = max(cl$times, na.rm = TRUE),
      freq_min = min(cl$freqs, na.rm = TRUE),
      freq_max = max(cl$freqs, na.rm = TRUE),
      CoMtime = cl$CoMtime,
      CoMfreq = cl$CoMfreq,
      peakZscore = cl$peakZscore,
      sumZscore = cl$sumZscore,
      peakpAdj = if (!is.null(cl$peakpAdj)) cl$peakpAdj else NA_real_,
      minpAdj = if (length(pAdj_valid) == 0) NA_real_ else min(pAdj_valid),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Tidy per-pixel cluster data
#'
#' Returns a long data.frame with one row per cluster pixel. Useful for plotting
#' or post-hoc summaries.
#'
#' @param clusters list of clusters as returned by \code{\link{detect_clusters}}.
#' @return data.frame with columns \code{cluster}, \code{time}, \code{freq},
#'   \code{Zscore}, and \code{pAdj} (if present).
#' @export
clusters_pixels <- function(clusters) {
  if (is.null(clusters) || length(clusters) == 0) {
    return(data.frame(
      cluster = integer(0),
      time = numeric(0),
      freq = numeric(0),
      Zscore = numeric(0),
      pAdj = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  dfs <- lapply(seq_along(clusters), function(i) {
    cl <- clusters[[i]]
    n <- length(cl$Zscores)
    pAdj <- if (!is.null(cl$pAdj)) cl$pAdj else rep(NA_real_, n)
    data.frame(
      cluster = i,
      time = cl$times,
      freq = cl$freqs,
      Zscore = cl$Zscores,
      pAdj = pAdj,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}
