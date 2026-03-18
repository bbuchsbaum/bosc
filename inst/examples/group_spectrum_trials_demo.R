#!/usr/bin/env Rscript

# Dense-sampling grouped-spectrum demo inspired by the trial-to-spectrum
# workflow used in Biba et al. (2026), while remaining fully synthetic and
# package-native.

repo_root <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(dirname(dirname(sub("^--file=", "", file_arg[1])))), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

if (!requireNamespace("bosc", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install 'bosc' or 'pkgload' to run this example.", call. = FALSE)
  }
  pkgload::load_all(repo_root(), quiet = TRUE)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root(), quiet = TRUE)
}

suppressPackageStartupMessages(library(bosc))

set.seed(42)

n_subject <- 16
n_bin <- 28
n_trial_per_bin <- 12
fs <- 30
target_freq <- 8
soa_bins <- round(seq(200, 1100, by = 1000 / 30), 2)

trial_dat <- expand.grid(
  subject = paste0("s", seq_len(n_subject)),
  soa_nominal = soa_bins,
  rep = seq_len(n_trial_per_bin),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
trial_dat$trial_index <- ave(trial_dat$rep, trial_dat$subject, FUN = seq_along)
trial_dat$task <- rep(c(-0.5, 0.5), length.out = nrow(trial_dat))
trial_dat$soa_measured <- trial_dat$soa_nominal + stats::runif(nrow(trial_dat), min = -12, max = 12)

bin_index <- match(trial_dat$soa_nominal, soa_bins)
subject_phase <- rep(seq(0, pi / 3, length.out = n_subject), each = n_bin * n_trial_per_bin)
subject_bias <- rep(stats::rnorm(n_subject, sd = 0.25), each = n_bin * n_trial_per_bin)
theta_signal <- sin(2 * pi * target_freq * (bin_index - 1) / fs + subject_phase)

trial_dat$value <- 2 +
  subject_bias +
  0.45 * trial_dat$task -
  0.01 * scale(trial_dat$trial_index)[, 1] +
  0.9 * theta_signal +
  stats::rnorm(nrow(trial_dat), sd = 0.35)

trial_dat <- align_time_bins(
  trial_dat,
  time = "soa_measured",
  bins = soa_bins,
  tolerance = 1000 / 60,
  out = "soa_bin"
)

trial_dat <- residualize_trials(
  trial_dat,
  response = "value",
  terms = c("task", "trial_index"),
  by = "subject",
  family = "gaussian",
  type = "response",
  out = "value_resid"
)

res <- group_spectrum_test_trials(
  trial_dat,
  value = "value_resid",
  bin = "soa_bin",
  subject = "subject",
  bins = soa_bins,
  fs = fs,
  flim = c(3, 14),
  null = "shuffle_bins",
  nrep = 200,
  detrend = TRUE,
  detrend_order = 2,
  pad_to = 30,
  taper = "hann",
  seed = 123
)

peak_idx <- which.max(res$observed$power)
peak_freq <- res$observed$freq[peak_idx]
peak_p <- res$stats$p_adj[peak_idx]

print(data.frame(
  freq = res$observed$freq,
  power = round(res$observed$power, 4),
  p_adj = signif(res$stats$p_adj, 3)
))

cat("\nPeak frequency:", peak_freq, "Hz\n")
cat("FDR-adjusted p at peak:", signif(peak_p, 3), "\n")

stopifnot(isTRUE(all.equal(res$observed$freq, 3:14)))
stopifnot(abs(peak_freq - target_freq) < 1e-12)
stopifnot(is.finite(peak_p), peak_p < 0.05)
