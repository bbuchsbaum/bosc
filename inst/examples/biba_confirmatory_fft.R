#!/usr/bin/env Rscript

# Reproduces the confirmatory stationary FFT results reported in:
# Biba, T. M., Decker, A., Herrmann, B., Fukuda, K., Katz, C. N.,
# Valiante, T. A., & Duncan, K. (2026). Episodic memory encoding
# fluctuates at a theta rhythm of 3-10 Hz. Nature Human Behaviour.
# https://doi.org/10.1038/s41562-026-02416-5

extra_libs <- Sys.getenv("BOSC_BIBA_R_LIBS", unset = "")
if (nzchar(extra_libs)) {
  .libPaths(c(strsplit(extra_libs, .Platform$path.sep)[[1]], .libPaths()))
}

required_pkgs <- c(
  "knitr", "dplyr", "readr", "ggplot2", "tibble", "magrittr",
  "stringr", "purrr", "tidyr", "DescTools", "prospectr", "e1071"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Install required packages before running this example: ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(magrittr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

message_line <- function(...) {
  cat(paste0(..., "\n"))
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE))
  }
  NULL
}

locate_repo_root <- function() {
  candidates <- c(
    getwd(),
    dirname(script_path() %||% getwd()),
    dirname(dirname(script_path() %||% getwd()))
  )
  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  hits <- candidates[
    file.exists(file.path(candidates, "DESCRIPTION")) &
      dir.exists(file.path(candidates, "MBO_shared", "Code"))
  ]
  if (length(hits) == 0) {
    stop(
      "Could not locate the bosc source root. Run this from the package root ",
      "or keep the script under inst/examples in the source checkout.",
      call. = FALSE
    )
  }
  hits[1]
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest)) {
    return(invisible(dest))
  }
  utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  invisible(dest)
}

load_bosc <- function(repo_root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_root, quiet = TRUE)
  } else if (!requireNamespace("bosc", quietly = TRUE)) {
    stop("Install 'bosc' or 'pkgload' to run this example.", call. = FALSE)
  }
  suppressPackageStartupMessages(library(bosc))
}

ensure_cleaning_helper <- function(code_dir, cache_dir) {
  clean_rmd <- file.path(code_dir, "MBO_data_cleaning_helper.Rmd")
  clean_r <- file.path(cache_dir, "clean_helper_exact.R")

  if (!file.exists(clean_r)) {
    knitr::purl(clean_rmd, output = clean_r, quiet = TRUE)
  }
  source(clean_r, local = .GlobalEnv)
}

assert_close <- function(value, target, tol, label) {
  if (!isTRUE(all.equal(unname(value), unname(target), tolerance = tol))) {
    stop(
      sprintf(
        "%s mismatch: got %.8f, expected %.8f (tol = %.8f)",
        label, as.numeric(value), as.numeric(target), tol
      ),
      call. = FALSE
    )
  }
}

extract_fft_row <- function(stats_df, hz) {
  row <- dplyr::filter(stats_df, freq == hz)
  if (nrow(row) != 1) {
    stop("Expected exactly one FFT row at ", hz, " Hz.", call. = FALSE)
  }
  row
}

stats_against_external_null <- function(observed_power, freqs, perm_path) {
  perm_df <- readr::read_csv(perm_path, show_col_types = FALSE)
  perm_freq <- suppressWarnings(as.numeric(names(perm_df)))
  if (anyNA(perm_freq)) {
    stop("Permutation CSV must use frequency-bin column names.", call. = FALSE)
  }
  idx <- match(freqs, perm_freq)
  if (anyNA(idx)) {
    stop("Observed frequencies did not match permutation CSV columns.", call. = FALSE)
  }
  perm_mat <- as.matrix(perm_df[, idx, drop = FALSE])
  z <- (observed_power - colMeans(perm_mat)) / apply(perm_mat, 2, stats::sd)
  p <- stats::pnorm(z, lower.tail = FALSE)
  tibble::tibble(
    freq = freqs,
    power = observed_power,
    z = z,
    p = p,
    p_adj = stats::p.adjust(p, method = "fdr")
  )
}

complete_subjects <- function(data, value, bins, extra_filter = NULL) {
  if (!is.null(extra_filter)) {
    data <- dplyr::filter(data, !!rlang::parse_expr(extra_filter))
  }
  data %>%
    dplyr::filter(!is.na(SOAcorr_ct), !is.na(.data[[value]])) %>%
    dplyr::count(subject, SOAcorr_ct) %>%
    dplyr::count(subject, name = "n_bin") %>%
    dplyr::filter(n_bin == length(bins)) %>%
    dplyr::pull(subject)
}

observed_group_spectrum <- function(data, value, bins, fun = mean) {
  series <- make_grouped_series(
    data = data,
    value = value,
    bin = "SOAcorr_ct",
    subject = "subject",
    bins = bins,
    fun = fun,
    incomplete = "drop"
  )
  series$data <- t(apply(series$data, 1, function(ts) {
    prospectr::detrend(ts, wav = seq_along(ts), p = 2) *
      e1071::hanning.window(length(ts))
  }))
  group_spectrum(
    series,
    fs = 30,
    flim = c(3, 14),
    detrend = FALSE,
    pad_to = 30,
    spectrum = "raw_power",
    taper = "none"
  )
}

observed_ies_spectrum <- function(data, hr_value, rt_value, rt_filter, bins) {
  hr_df <- aggregate_by_bin(
    data = dplyr::filter(data, latency > 300),
    value = hr_value,
    bin = "SOAcorr_ct",
    by = "subject",
    bins = bins,
    fun = mean,
    na_rm = TRUE,
    complete = TRUE,
    out = "hr"
  ) %>%
    dplyr::select(subject, SOAcorr_ct, hr)

  rt_df <- aggregate_by_bin(
    data = dplyr::filter(data, !!rlang::parse_expr(rt_filter)),
    value = rt_value,
    bin = "SOAcorr_ct",
    by = "subject",
    bins = bins,
    fun = stats::median,
    na_rm = TRUE,
    complete = TRUE,
    out = "rt"
  ) %>%
    dplyr::select(subject, SOAcorr_ct, rt)

  ies_df <- dplyr::left_join(hr_df, rt_df, by = c("subject", "SOAcorr_ct")) %>%
    dplyr::mutate(ies = rt / hr)

  complete_subs <- ies_df %>%
    dplyr::group_by(subject) %>%
    dplyr::summarize(ok = n() == length(bins) && all(is.finite(ies)), .groups = "drop") %>%
    dplyr::filter(ok) %>%
    dplyr::pull(subject)

  ies_df <- dplyr::filter(ies_df, subject %in% complete_subs)
  subjects <- sort(unique(ies_df$subject))
  ies_mat <- t(vapply(subjects, function(subj) {
    dplyr::filter(ies_df, subject == subj) %>%
      dplyr::arrange(SOAcorr_ct) %>%
      dplyr::pull(ies)
  }, numeric(length(bins))))

  ies_series <- list(
    data = ies_mat,
    subjects = subjects,
    bins = bins,
    measure = "ies",
    meta = list(input = "derived_ratio")
  )
  ies_series$data <- t(apply(ies_series$data, 1, function(ts) {
    prospectr::detrend(ts, wav = seq_along(ts), p = 2) *
      e1071::hanning.window(length(ts))
  }))
  group_spectrum(
    ies_series,
    fs = 30,
    flim = c(3, 14),
    detrend = FALSE,
    pad_to = 30,
    spectrum = "raw_power",
    taper = "none"
  )
}

run_confirmatory_fft <- function(mem_resid, cls_resid, perm_paths) {
  bins <- sort(unique(mem_resid$SOAcorr_ct[!is.na(mem_resid$SOAcorr_ct)]))

  mem_hr_obs <- observed_group_spectrum(
    dplyr::filter(mem_resid, item_cor == "Old", latency > 300),
    value = "hr_resid_tr",
    bins = bins
  )
  mem_rt_obs <- observed_group_spectrum(
    dplyr::filter(mem_resid, asso_hit == 1, latency > 300),
    value = "rt_resid_tr",
    bins = bins,
    fun = stats::median
  )
  mem_ies_obs <- observed_ies_spectrum(
    dplyr::filter(mem_resid, item_cor == "Old"),
    hr_value = "hr_resid_tr",
    rt_value = "rt_resid_tr",
    rt_filter = "asso_hit == 1 & latency > 300",
    bins = bins
  )

  cls_hr_obs <- observed_group_spectrum(
    dplyr::filter(cls_resid, latency > 300),
    value = "hr_resid_tr",
    bins = bins
  )
  cls_rt_obs <- observed_group_spectrum(
    dplyr::filter(cls_resid, class_TS_hit == 1, latency > 300),
    value = "rt_resid_tr",
    bins = bins,
    fun = stats::median
  )
  cls_ies_obs <- observed_ies_spectrum(
    cls_resid,
    hr_value = "hr_resid_tr",
    rt_value = "rt_resid_tr",
    rt_filter = "class_TS_hit == 1 & latency > 300",
    bins = bins
  )

  list(
    mem_hr = stats_against_external_null(mem_hr_obs$observed$power, mem_hr_obs$observed$freq, perm_paths$mem_hr),
    mem_rt = stats_against_external_null(mem_rt_obs$observed$power, mem_rt_obs$observed$freq, perm_paths$mem_rt),
    mem_ies = stats_against_external_null(mem_ies_obs$observed$power, mem_ies_obs$observed$freq, perm_paths$mem_ies),
    cls_hr = stats_against_external_null(cls_hr_obs$observed$power, cls_hr_obs$observed$freq, perm_paths$cls_hr),
    cls_rt = stats_against_external_null(cls_rt_obs$observed$power, cls_rt_obs$observed$freq, perm_paths$cls_rt),
    cls_ies = stats_against_external_null(cls_ies_obs$observed$power, cls_ies_obs$observed$freq, perm_paths$cls_ies)
  )
}

repo_root <- locate_repo_root()
load_bosc(repo_root)

cache_dir <- Sys.getenv(
  "BOSC_BIBA_CACHE",
  unset = file.path(tempdir(), "bosc_biba_confirmatory")
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

urls <- list(
  raw = "https://osf.io/download/k5zxp/",
  survey = "https://osf.io/download/t8k4v/",
  mem_hr = "https://osf.io/download/6834d3e85cb203bc0f85b52f/",
  mem_rt = "https://osf.io/download/6834d3eb32ad3ecae9a4dd1b/",
  mem_ies = "https://osf.io/download/6834d3e98f795288f2ab8511/",
  cls_hr = "https://osf.io/download/6834d3e43c1fd3afcda4dc03/",
  cls_rt = "https://osf.io/download/6834d3e78f795288f2ab850f/",
  cls_ies = "https://osf.io/download/6834d3e41636412fbdd6a499/"
)

paths <- list(
  raw = file.path(cache_dir, "confirmatory_raw.csv"),
  survey = file.path(cache_dir, "confirmatory_survey.csv"),
  mem_hr = file.path(cache_dir, "FFTsur_mem_hr.csv"),
  mem_rt = file.path(cache_dir, "FFTsur_mem_rt.csv"),
  mem_ies = file.path(cache_dir, "FFTsur_mem_ies.csv"),
  cls_hr = file.path(cache_dir, "FFTsur_cls_hr.csv"),
  cls_rt = file.path(cache_dir, "FFTsur_cls_rt.csv"),
  cls_ies = file.path(cache_dir, "FFTsur_cls_ies.csv")
)

message_line("Using bosc root: ", repo_root)
message_line("Using cache dir: ", cache_dir)
purrr::walk2(urls, paths, download_if_missing)

ensure_cleaning_helper(file.path(repo_root, "MBO_shared", "Code"), cache_dir)

# The confirmatory FFT does not use questionnaire fields, but the original
# helper assumes every subject has a survey row. The OSF export does not.
bindQuesInf <- function(inEnc, inRet, inQues, colsurg) {
  subject_id <- unique(inEnc$subject)
  ques <- dplyr::filter(inQues, subject == subject_id)

  if (nrow(ques) == 0) {
    enc_ques <- as.data.frame(matrix(NA, nrow = nrow(inEnc), ncol = length(colsurg)))
    ret_ques <- as.data.frame(matrix(NA, nrow = nrow(inRet), ncol = length(colsurg)))
    colnames(enc_ques) <- colsurg
    colnames(ret_ques) <- colsurg
  } else {
    ques <- ques[1, colsurg, drop = FALSE]
    enc_ques <- ques[rep(1, nrow(inEnc)), , drop = FALSE]
    ret_ques <- ques[rep(1, nrow(inRet)), , drop = FALSE]
  }

  list(cbind(inEnc, enc_ques), cbind(inRet, ret_ques))
}

resid_rds <- file.path(cache_dir, "confirmatory_residualized_exact.rds")
if (file.exists(resid_rds)) {
  resid_obj <- readRDS(resid_rds)
  if (!all(c("cls_resid", "mem_resid") %in% names(resid_obj))) {
    if (all(c("enc", "ret") %in% names(resid_obj))) {
      if (all(c("hr_resid_tr", "rt_resid_tr") %in% names(resid_obj$enc)) &&
          all(c("hr_resid_tr", "rt_resid_tr") %in% names(resid_obj$ret))) {
        resid_obj <- list(
          cls_resid = resid_obj$enc,
          mem_resid = resid_obj$ret
        )
      } else {
        resid_obj <- getResid_ind(
          inEnc = resid_obj$enc,
          inRet = resid_obj$ret,
          rmTri1 = TRUE,
          wPlots = FALSE,
          lat = 300
        )
      }
      saveRDS(resid_obj, resid_rds)
    } else {
      stop("Existing residual cache has an unrecognized structure.", call. = FALSE)
    }
  }
} else {
  raw_dat <- as.data.frame(readr::read_csv(paths$raw, show_col_types = FALSE))
  survey_dat <- as.data.frame(readr::read_csv(paths$survey, show_col_types = FALSE))
  stim_list <- as.data.frame(readr::read_csv(
    file.path(repo_root, "MBO_shared", "Stimuli", "full_dt_wdr.csv"),
    show_col_types = FALSE
  ))

  cleaned <- datCLWrap(raw_dat, survey_dat, stim_list)
  perf <- chancePerf(cleaned, plots = FALSE)
  keep_subjects <- perf[[2]]$subject

  enc_dat <- dplyr::filter(cleaned[[1]], subject %in% keep_subjects)
  ret_dat <- dplyr::filter(cleaned[[2]], subject %in% keep_subjects)

  resid_obj <- getResid_ind(
    inEnc = enc_dat,
    inRet = ret_dat,
    rmTri1 = TRUE,
    wPlots = FALSE,
    lat = 300
  )
  saveRDS(resid_obj, resid_rds)
}

cls_resid <- resid_obj$cls_resid
mem_resid <- resid_obj$mem_resid
n_subject <- length(unique(mem_resid$subject))
if (!identical(n_subject, 125L)) {
  stop("Expected 125 confirmatory subjects after exclusions, got ", n_subject, ".", call. = FALSE)
}

fft_results <- run_confirmatory_fft(mem_resid, cls_resid, paths[c(
  "mem_hr", "mem_rt", "mem_ies", "cls_hr", "cls_rt", "cls_ies"
)])

mem_hr_7 <- extract_fft_row(fft_results$mem_hr, 7)
mem_rt_13 <- extract_fft_row(fft_results$mem_rt, 13)
mem_ies_7 <- extract_fft_row(fft_results$mem_ies, 7)

assert_close(mem_hr_7$z, 2.644, 0.001, "Memory HR z at 7 Hz")
assert_close(mem_hr_7$p, 0.004095484, 1e-6, "Memory HR raw p at 7 Hz")
assert_close(mem_hr_7$p_adj, 0.0491458, 1e-6, "Memory HR FDR p at 7 Hz")

assert_close(mem_ies_7$z, 2.688, 0.001, "Memory IES z at 7 Hz")
assert_close(mem_ies_7$p, 0.003595549, 1e-6, "Memory IES raw p at 7 Hz")
assert_close(mem_ies_7$p_adj, 0.04314659, 1e-6, "Memory IES FDR p at 7 Hz")

assert_close(mem_rt_13$p, 0.05966694, 1e-6, "Memory RT raw p at 13 Hz")

cls_min_fdr <- c(
  min(fft_results$cls_hr$p_adj, na.rm = TRUE),
  min(fft_results$cls_rt$p_adj, na.rm = TRUE),
  min(fft_results$cls_ies$p_adj, na.rm = TRUE)
)
if (any(cls_min_fdr <= 0.05)) {
  stop("Classification control FFT unexpectedly became significant after FDR.", call. = FALSE)
}

summary_tbl <- tibble::tibble(
  measure = c("mem_hr", "mem_rt", "mem_ies", "cls_hr", "cls_rt", "cls_ies"),
  peak_hz = c(
    fft_results$mem_hr$freq[which.min(fft_results$mem_hr$p)],
    fft_results$mem_rt$freq[which.min(fft_results$mem_rt$p)],
    fft_results$mem_ies$freq[which.min(fft_results$mem_ies$p)],
    fft_results$cls_hr$freq[which.min(fft_results$cls_hr$p)],
    fft_results$cls_rt$freq[which.min(fft_results$cls_rt$p)],
    fft_results$cls_ies$freq[which.min(fft_results$cls_ies$p)]
  ),
  min_raw_p = c(
    min(fft_results$mem_hr$p, na.rm = TRUE),
    min(fft_results$mem_rt$p, na.rm = TRUE),
    min(fft_results$mem_ies$p, na.rm = TRUE),
    min(fft_results$cls_hr$p, na.rm = TRUE),
    min(fft_results$cls_rt$p, na.rm = TRUE),
    min(fft_results$cls_ies$p, na.rm = TRUE)
  ),
  min_fdr_p = c(
    min(fft_results$mem_hr$p_adj, na.rm = TRUE),
    min(fft_results$mem_rt$p_adj, na.rm = TRUE),
    min(fft_results$mem_ies$p_adj, na.rm = TRUE),
    min(fft_results$cls_hr$p_adj, na.rm = TRUE),
    min(fft_results$cls_rt$p_adj, na.rm = TRUE),
    min(fft_results$cls_ies$p_adj, na.rm = TRUE)
  )
)

out_csv <- file.path(cache_dir, "biba_confirmatory_fft_summary.csv")
readr::write_csv(summary_tbl, out_csv)

message_line("")
message_line("Confirmatory sample size: ", n_subject)
print(summary_tbl)
message_line("")
message_line("Wrote summary CSV: ", out_csv)
message_line("All confirmatory checks passed.")
