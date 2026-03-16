#' Align observed times to a target bin grid
#'
#' Maps observed timing values to the nearest target bin. This is useful when
#' nominal bins are known in advance but measured times contain small timing
#' deviations.
#'
#' @param data data.frame containing the timing column.
#' @param time character scalar naming the observed time column.
#' @param bins numeric vector of target bin values.
#' @param tolerance optional non-negative scalar. If provided, values farther
#'   than \code{tolerance} from the nearest target bin are assigned \code{NA}.
#' @param method matching method. Currently only \code{"nearest"} is supported.
#' @param out optional output column name. Defaults to \code{paste0(time, "_bin")}.
#' @return data.frame with an added aligned-bin column.
#' @export
align_time_bins <- function(data,
                            time,
                            bins,
                            tolerance = NULL,
                            method = c("nearest"),
                            out = NULL) {
  data <- validate_data_frame(data)
  time <- validate_column_name(time, "time")
  if (!time %in% names(data)) stop("time column not found in data.")
  bins <- validate_bins(bins)
  method <- match.arg(method)
  if (!is.null(tolerance)) {
    if (length(tolerance) != 1 || !is.finite(tolerance) || tolerance < 0) {
      stop("tolerance must be NULL or a non-negative scalar.")
    }
  }
  if (is.null(out)) out <- paste0(time, "_bin")
  out <- validate_column_name(out, "out")
  if (method != "nearest") stop("Unsupported method.")

  vals <- data[[time]]
  aligned <- rep(NA_real_, length(vals))
  ok <- is.finite(vals)
  if (any(ok)) {
    nearest_idx <- vapply(vals[ok], function(x) {
      which.min(abs(bins - x))
    }, integer(1))
    nearest <- bins[nearest_idx]
    if (!is.null(tolerance)) {
      keep <- abs(vals[ok] - nearest) <= tolerance
      nearest[!keep] <- NA_real_
    }
    aligned[ok] <- nearest
  }

  data[[out]] <- aligned
  data
}

#' Residualize trial-wise responses against nuisance terms
#'
#' Fits linear models within optional groups and returns response residuals in a
#' new column. This provides a general trial-level nuisance-regression step for
#' behavioral analyses without hard-coding any specific covariates.
#'
#' @param data data.frame containing response and nuisance columns.
#' @param response character scalar naming the response column.
#' @param terms optional character vector of nuisance column names. If empty or
#'   \code{NULL}, the response is copied unchanged.
#' @param by optional character vector of grouping columns. Separate models are
#'   fit within each group.
#' @param family model family. May be either a character scalar such as
#'   \code{"gaussian"}, \code{"binomial"}, or \code{"inverse.gaussian"}, or a
#'   family object returned by \code{\link[stats]{family}}.
#' @param type residual type passed to \code{\link[stats]{residuals.glm}}.
#' @param out optional output column name. Defaults to \code{paste0(response,
#'   "_resid")}.
#' @return data.frame with an added residual column.
#' @export
residualize_trials <- function(data,
                               response,
                               terms = NULL,
                               by = NULL,
                               family = "gaussian",
                               type = c("response", "working", "deviance", "pearson"),
                               out = NULL) {
  data <- validate_data_frame(data)
  response <- validate_column_name(response, "response")
  if (!response %in% names(data)) stop("response column not found in data.")
  terms <- validate_optional_columns(terms, data, "terms")
  by <- validate_optional_columns(by, data, "by")
  type <- match.arg(type)
  family_obj <- validate_model_family(family)
  if (is.null(out)) out <- paste0(response, "_resid")
  out <- validate_column_name(out, "out")

  resid_all <- rep(NA_real_, nrow(data))
  if (length(terms) == 0) {
    resid_all <- as.numeric(data[[response]])
    data[[out]] <- resid_all
    return(data)
  }

  group_index <- split_row_index(data, by)
  form <- stats::as.formula(
    paste(response, "~", paste(terms, collapse = " + "))
  )

  for (idx in group_index) {
    df <- data[idx, c(response, terms), drop = FALSE]
    cc <- stats::complete.cases(df)
    if (!any(cc)) next
    fit <- try(
      stats::glm(
        form,
        data = df[cc, , drop = FALSE],
        family = family_obj,
        na.action = stats::na.exclude
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    resid_all[idx[cc]] <- as.numeric(stats::residuals(fit, type = type))
  }

  data[[out]] <- resid_all
  data
}

#' Build a grouped-series object from trial-level data
#'
#' Aggregates trial values into regular bins for each subject or
#' subject-by-condition group, then reshapes the result into the grouped-series
#' structure accepted by \code{\link{group_spectrum}} and related functions.
#'
#' @param data data.frame containing the trial-wise values.
#' @param value character scalar naming the trial-level value column.
#' @param bin character scalar naming the bin column.
#' @param subject character scalar naming the subject column.
#' @param by optional character vector of additional grouping columns.
#' @param bins optional vector of target bins.
#' @param fun summary function applied within each bin.
#' @param na_rm logical; if \code{TRUE}, remove missing trial values before
#'   aggregation.
#' @param complete logical; if \code{TRUE}, request explicit rows for missing
#'   group/bin combinations.
#' @param incomplete handling for subject groups with missing or non-finite bin
#'   summaries. \code{"error"} stops, while \code{"drop"} removes incomplete
#'   groups.
#' @param out optional output measure name. Defaults to \code{value}.
#' @param meta optional named list of metadata to attach.
#' @return grouped-series list with \code{data}, \code{subjects}, \code{bins},
#'   \code{measure}, \code{meta}, plus \code{aggregated} and \code{group_keys}.
#' @keywords internal
make_grouped_series <- function(data,
                                value,
                                bin,
                                subject,
                                by = NULL,
                                bins = NULL,
                                fun = mean,
                                na_rm = TRUE,
                                complete = TRUE,
                                incomplete = c("error", "drop"),
                                out = NULL,
                                meta = list()) {
  data <- validate_data_frame(data)
  value <- validate_column_name(value, "value")
  bin <- validate_column_name(bin, "bin")
  subject <- validate_column_name(subject, "subject")
  if (!value %in% names(data)) stop("value column not found in data.")
  if (!bin %in% names(data)) stop("bin column not found in data.")
  if (!subject %in% names(data)) stop("subject column not found in data.")
  by <- validate_optional_columns(by, data, "by")
  incomplete <- match.arg(incomplete)
  if (!is.list(meta)) stop("meta must be a list.")
  group_cols <- c(subject, by)

  agg <- aggregate_by_bin(
    data = data,
    value = value,
    bin = bin,
    by = group_cols,
    bins = bins,
    fun = fun,
    na_rm = na_rm,
    complete = complete,
    out = out
  )
  if (nrow(agg) == 0) {
    stop("No grouped data available after aggregation.")
  }

  if (is.null(bins)) {
    bins_final <- sort(unique(agg[[bin]][!is.na(agg[[bin]])]))
  } else {
    bins_final <- sort(unique(bins))
  }
  measure_name <- out %||% value
  if (length(bins_final) < 2) {
    stop("At least two bins are required to form a grouped series.")
  }

  key_df <- unique_group_frame(agg, group_cols)
  is_complete <- vapply(seq_len(nrow(key_df)), function(i) {
    rows <- subset_group_rows(agg, key_df, i, group_cols)
    idx <- match(bins_final, rows[[bin]])
    if (any(is.na(idx))) return(FALSE)
    vals <- rows[[measure_name]][idx]
    all(is.finite(vals))
  }, logical(1))

  if (!all(is_complete)) {
    if (incomplete == "error") {
      bad <- group_labels(key_df[!is_complete, , drop = FALSE], group_cols)
      stop(
        "Incomplete grouped series for: ",
        paste(bad, collapse = ", ")
      )
    }
    key_df <- key_df[is_complete, , drop = FALSE]
  }
  if (nrow(key_df) == 0) {
    stop("No complete subject groups remain after filtering.")
  }

  mat <- matrix(NA_real_, nrow = nrow(key_df), ncol = length(bins_final))
  for (i in seq_len(nrow(key_df))) {
    rows <- subset_group_rows(agg, key_df, i, group_cols)
    idx <- match(bins_final, rows[[bin]])
    mat[i, ] <- rows[[measure_name]][idx]
  }
  rownames(mat) <- group_labels(key_df, group_cols)

  list(
    data = mat,
    subjects = rownames(mat),
    bins = bins_final,
    measure = measure_name,
    meta = c(
      list(
        input = "trial_data",
        subject = subject,
        by = by,
        bin = bin,
        incomplete = incomplete
      ),
      meta
    ),
    aggregated = agg,
    group_keys = key_df
  )
}

#' Permute trial labels within groups
#'
#' Shuffles or circularly shifts a label column within groups. This is useful
#' for trial-level permutation nulls where the trial values are left unchanged
#' but their bin assignments are resampled.
#'
#' @param data data.frame containing the label column.
#' @param column character scalar naming the column to permute.
#' @param by optional character vector of grouping columns.
#' @param method permutation method. \code{"shuffle"} samples without
#'   replacement. \code{"circular_shift"} preserves local order after sorting by
#'   \code{order_by}.
#' @param order_by optional character scalar naming the within-group ordering
#'   column required for \code{"circular_shift"}.
#' @param out optional output column name. Defaults to \code{column}.
#' @param seed optional random seed.
#' @return data.frame with the permuted label column.
#' @keywords internal
permute_trial_labels <- function(data,
                                 column,
                                 by = NULL,
                                 method = c("shuffle", "circular_shift"),
                                 order_by = NULL,
                                 out = NULL,
                                 seed = NULL) {
  data <- validate_data_frame(data)
  column <- validate_column_name(column, "column")
  if (!column %in% names(data)) stop("column not found in data.")
  by <- validate_optional_columns(by, data, "by")
  method <- match.arg(method)
  if (is.null(out)) out <- column
  out <- validate_column_name(out, "out")
  if (method == "circular_shift") {
    order_by <- validate_column_name(order_by, "order_by")
    if (!order_by %in% names(data)) stop("order_by column not found in data.")
  }

  permuted <- data[[column]]
  group_index <- split_row_index(data, by)

  permuted <- with_seed(seed, {
    out_vals <- permuted
    for (idx in group_index) {
      if (length(idx) <= 1L) next
      vals <- data[[column]][idx]
      if (method == "shuffle") {
        out_vals[idx] <- sample(vals, size = length(vals), replace = FALSE)
      } else {
        ord <- order(data[[order_by]][idx], seq_along(idx), na.last = TRUE)
        shift <- sample.int(length(idx), size = 1L) - 1L
        out_vals[idx[ord]] <- circ_shift_vector(vals[ord], shift)
      }
    }
    out_vals
  })

  data[[out]] <- permuted
  data
}

#' Aggregate trial values into regularized bins
#'
#' Aggregates trial-level values within bins and optional grouping variables,
#' producing a tidy table suitable for grouped time-series analysis.
#'
#' @param data data.frame containing the value and bin columns.
#' @param value character scalar naming the value column to aggregate.
#' @param bin character scalar naming the bin column.
#' @param by optional character vector of grouping columns.
#' @param bins optional vector of target bins. If \code{NULL}, uses the sorted
#'   unique non-missing values observed in \code{data[[bin]]}.
#' @param fun summary function applied within each bin.
#' @param na_rm logical; if \code{TRUE}, remove missing values before applying
#'   \code{fun}.
#' @param complete logical; if \code{TRUE}, return explicit rows for missing
#'   group/bin combinations with \code{NA} summaries.
#' @param out optional output column name. Defaults to \code{value}.
#' @return data.frame with grouping columns, \code{bin}, the aggregated value,
#'   and count columns \code{n} and \code{n_nonmissing}.
#' @export
aggregate_by_bin <- function(data,
                             value,
                             bin,
                             by = NULL,
                             bins = NULL,
                             fun = mean,
                             na_rm = TRUE,
                             complete = TRUE,
                             out = NULL) {
  data <- validate_data_frame(data)
  value <- validate_column_name(value, "value")
  bin <- validate_column_name(bin, "bin")
  if (!value %in% names(data)) stop("value column not found in data.")
  if (!bin %in% names(data)) stop("bin column not found in data.")
  by <- validate_optional_columns(by, data, "by")
  if (is.null(out)) out <- value
  out <- validate_column_name(out, "out")
  if (!is.function(fun)) stop("fun must be a function.")

  if (is.null(bins)) {
    bins <- sort(unique(data[[bin]][!is.na(data[[bin]])]))
  } else {
    bins <- sort(unique(bins))
  }

  if (length(bins) == 0) {
    out_df <- unique_group_frame(data, by)
    if (!complete || nrow(out_df) == 0) {
      return(data.frame())
    }
    out_df[[bin]] <- vector(mode = typeof(data[[bin]]), length = 0)
    out_df[[out]] <- numeric(0)
    out_df$n <- integer(0)
    out_df$n_nonmissing <- integer(0)
    return(out_df)
  }

  groups <- unique_group_frame(data, by)
  if (nrow(groups) == 0) groups <- data.frame(.row = 1L)[, FALSE, drop = FALSE]

  rows <- vector("list", nrow(groups) * if (complete) length(bins) else max(1L, length(bins)))
  row_i <- 0L

  for (g in seq_len(nrow(groups))) {
    group_mask <- rep(TRUE, nrow(data))
    if (length(by) > 0) {
      for (col in by) {
        group_mask <- group_mask & same_or_na(data[[col]], groups[[col]][g])
      }
    }

    group_data <- data[group_mask, , drop = FALSE]
    group_bins <- if (complete) bins else sort(unique(group_data[[bin]][!is.na(group_data[[bin]])]))
    for (b in group_bins) {
      in_bin <- !is.na(group_data[[bin]]) & group_data[[bin]] == b
      vals <- group_data[[value]][in_bin]
      n_total <- length(vals)
      n_nonmissing <- sum(!is.na(vals))
      if (na_rm) vals <- vals[!is.na(vals)]
      agg <- if (length(vals) == 0) NA_real_ else as.numeric(fun(vals))
      row_i <- row_i + 1L
      row <- groups[g, , drop = FALSE]
      row[[bin]] <- b
      row[[out]] <- agg
      row$n <- n_total
      row$n_nonmissing <- n_nonmissing
      rows[[row_i]] <- row
    }
  }

  out_df <- do.call(rbind, rows[seq_len(row_i)])
  rownames(out_df) <- NULL
  out_df
}

validate_data_frame <- function(data) {
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  data
}

validate_column_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1 || !nzchar(x)) {
    stop(arg, " must be a non-empty character scalar.")
  }
  x
}

validate_optional_columns <- function(x, data, arg) {
  if (is.null(x)) return(character(0))
  if (!is.character(x)) stop(arg, " must be NULL or a character vector.")
  if (length(x) == 0) return(character(0))
  missing_cols <- setdiff(x, names(data))
  if (length(missing_cols) > 0) {
    stop(arg, " columns not found in data: ", paste(missing_cols, collapse = ", "))
  }
  x
}

validate_bins <- function(bins) {
  if (!is.numeric(bins) || length(bins) == 0) {
    stop("bins must be a non-empty numeric vector.")
  }
  bins <- sort(unique(as.numeric(bins)))
  bins[is.finite(bins)]
}

validate_model_family <- function(family) {
  if (is.character(family) && length(family) == 1L) {
    fam_name <- family
    family <- switch(
      fam_name,
      gaussian = stats::gaussian(),
      binomial = stats::binomial(),
      poisson = stats::poisson(),
      Gamma = stats::Gamma(),
      inverse.gaussian = stats::inverse.gaussian(),
      stop("Unsupported family: ", fam_name)
    )
  }
  if (!is.list(family) || is.null(family$family) || !is.function(family$linkfun)) {
    stop("family must be a supported character scalar or a stats family object.")
  }
  family
}

split_row_index <- function(data, by) {
  if (length(by) == 0) return(list(seq_len(nrow(data))))
  grp <- interaction(data[by], drop = TRUE, lex.order = TRUE)
  split(seq_len(nrow(data)), grp)
}

unique_group_frame <- function(data, by) {
  if (length(by) == 0) return(data.frame(.group = 1L)[, FALSE, drop = FALSE])
  out <- unique(data[by], incomparables = FALSE)
  ord <- do.call(order, unname(out))
  out[ord, , drop = FALSE]
}

same_or_na <- function(x, y) {
  if (is.na(y)) {
    is.na(x)
  } else {
    x == y
  }
}

subset_group_rows <- function(data, keys, i, cols) {
  mask <- rep(TRUE, nrow(data))
  if (length(cols) > 0) {
    for (col in cols) {
      mask <- mask & same_or_na(data[[col]], keys[[col]][i])
    }
  }
  data[mask, , drop = FALSE]
}

group_labels <- function(keys, cols) {
  if (length(cols) == 1L) return(as.character(keys[[cols]]))
  apply(keys[cols], 1, function(x) paste(as.character(x), collapse = ":"))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
