make_group_sine <- function(n_subject = 4, n_bin = 80, fs = 24, freq = 6) {
  t <- seq(0, by = 1 / fs, length.out = n_bin)
  x <- vapply(seq_len(n_subject), function(i) {
    sin(2 * pi * freq * t + (i - 1) * pi / 5)
  }, numeric(n_bin))
  t(x)
}

test_that("as_grouped_series validates matrix and list inputs", {
  as_gs <- getFromNamespace("as_grouped_series", "bosc")
  expect_error(as_gs(matrix("a", 2, 3)), "numeric")
  expect_error(as_gs(matrix(1:2, nrow = 2)), "two bins")
  expect_error(as_gs(list(data = matrix(1:2, nrow = 2))), "two bins")
  expect_error(as_gs(list(data = matrix(1:6, 2), subjects = c("a"))), "subjects")
  expect_error(as_gs(list(data = matrix(1:6, 2), bins = 1:2)), "bins")
  expect_error(as_gs(list(x = 1:6)), "numeric matrix or grouped-series")
})

test_that("group_spectrum honors average=FALSE and validates sampling args", {
  x <- make_group_sine()
  res <- group_spectrum(x, fs = 24, flim = c(3, 10), average = FALSE, taper = "hann")
  expect_null(res$observed$power)
  expect_true(is.matrix(res$observed$per_subject))

  expect_error(group_spectrum(x, fs = 24, flim = c(10, 3)), "increasing positive")
  expect_error(group_spectrum(x, fs = 24, detrend_order = -1), "non-negative")
  expect_error(group_spectrum_test(x, fs = 24, nrep = 5, seed = c(1, 2)), "seed")
})

test_that("group_spectrum_test handles single null replicate", {
  x <- make_group_sine(n_subject = 3, n_bin = 60)
  res <- group_spectrum_test(
    x,
    fs = 24,
    flim = c(3, 10),
    null = "shuffle_labels",
    nrep = 1,
    seed = 11
  )
  expect_equal(dim(res$null$power), c(1L, length(res$observed$freq)))
  expect_true(all(is.finite(res$stats$p)))
})

test_that("group_spectrum_test_trials validates trial table inputs", {
  gst <- getFromNamespace("group_spectrum_test_trials", "bosc")
  dat <- data.frame(subject = "s1", bin = 1:3, value = 1:3)
  expect_error(
    gst(dat, value = "missing", bin = "bin", subject = "subject", fs = 10),
    "value column"
  )
  expect_error(
    gst(dat, value = "value", bin = "missing", subject = "subject", fs = 10),
    "bin column"
  )
  expect_error(
    gst(
      dat,
      value = "value",
      bin = "bin",
      subject = "subject",
      fs = 10,
      null = "circular_shift_bins"
    ),
    "order_by"
  )
  expect_error(
    gst(
      dat,
      value = "value",
      bin = "bin",
      subject = "subject",
      fs = 10,
      max_resample = 0
    ),
    "max_resample"
  )
})

test_that("group_spectrum_test_trials supports circular shift and resample repair", {
  gst <- getFromNamespace("group_spectrum_test_trials", "bosc")
  set.seed(301)
  x <- make_group_sine(n_subject = 3, n_bin = 48, fs = 24, freq = 6)
  trial_dat <- data.frame(
    subject = rep(paste0("s", seq_len(nrow(x))), each = ncol(x)),
    trial = rep(seq_len(ncol(x)), times = nrow(x)),
    bin = rep(seq_len(ncol(x)), times = nrow(x)),
    value = as.numeric(t(x))
  )

  circ <- gst(
    trial_dat,
    value = "value",
    bin = "bin",
    subject = "subject",
    order_by = "trial",
    fs = 24,
    bins = seq_len(ncol(x)),
    flim = c(3, 10),
    null = "circular_shift_bins",
    nrep = 4,
    seed = 88
  )
  expect_equal(circ$null$method, "circular_shift_bins")
  expect_true(all(is.finite(circ$stats$p)))

  sparse <- trial_dat
  sparse$value[sparse$subject == "s1" & sparse$bin == 2] <- NA_real_
  resampled <- gst(
    sparse,
    value = "value",
    bin = "bin",
    subject = "subject",
    fs = 24,
    bins = seq_len(ncol(x)),
    flim = c(3, 10),
    null = "shuffle_bins",
    incomplete = "drop",
    repair_incomplete = "resample",
    nrep = 3,
    max_resample = 50,
    seed = 99
  )
  expect_equal(resampled$meta$repair_incomplete, "resample")
  expect_true(all(is.finite(resampled$stats$p)))
})

test_that("internal grouped helpers cover edge branches", {
  circ_shift <- getFromNamespace("circ_shift_vector", "bosc")
  expect_equal(circ_shift(numeric(0), 1), numeric(0))
  expect_equal(circ_shift(c(1, 2, 3), 0), c(1, 2, 3))
  expect_equal(circ_shift(c(1, 2, 3), 3), c(1, 2, 3))

  with_seed <- getFromNamespace("with_seed", "bosc")
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  draw1 <- with_seed(123, rnorm(3))
  draw2 <- with_seed(123, rnorm(3))
  expect_equal(draw1, draw2)
  if (old_exists) {
    expect_equal(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), old_seed)
  } else {
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  }

  filter_keys <- getFromNamespace("filter_to_group_keys", "bosc")
  dat <- data.frame(subject = c("s1", "s1", "s2"), cond = c("a", "a", "b"), value = 1:3)
  keys <- data.frame(subject = c("s1", "s2"), cond = c("a", "b"))
  filtered <- filter_keys(dat, keys, c("subject", "cond"))
  expect_equal(nrow(filtered), 3L)
})

test_that("make_grouped_series and aggregate_by_bin validate edge cases", {
  mgs <- getFromNamespace("make_grouped_series", "bosc")
  dat <- data.frame(subject = "s1", bin = 1, value = 1)
  expect_error(
    mgs(dat, value = "value", bin = "bin", subject = "subject", bins = 1),
    "At least two bins"
  )

  incomplete <- data.frame(
    subject = c("s1", "s1", "s2", "s2"),
    bin = c(1, 2, 1, 3),
    value = c(10, 20, 11, 31)
  )
  expect_error(
    mgs(
      incomplete,
      value = "value",
      bin = "bin",
      subject = "subject",
      bins = 1:3,
      incomplete = "error"
    ),
    "Incomplete grouped series"
  )

  expect_error(aggregate_by_bin(dat, value = "value", bin = "bin", fun = "mean"), "function")
  expect_error(getFromNamespace("validate_data_frame", "bosc")(list(a = 1)), "data.frame")
  expect_error(getFromNamespace("validate_column_name", "bosc")("", "time"), "non-empty")
  expect_error(getFromNamespace("validate_optional_columns", "bosc")("missing", dat, "by"), "not found")
  expect_error(getFromNamespace("validate_bins", "bosc")(numeric(0)), "non-empty")
  expect_error(getFromNamespace("validate_model_family", "bosc")("badfamily"), "Unsupported family")
})

test_that("align_time_bins and permute_trial_labels validate arguments", {
  ptl <- getFromNamespace("permute_trial_labels", "bosc")
  dat <- data.frame(time = 1:3)
  expect_error(align_time_bins(dat, time = "time", bins = 1:3, tolerance = -1), "tolerance")
  expect_error(
    align_time_bins(dat, time = "time", bins = 1:3, method = "nearest", out = ""),
    "non-empty"
  )

  trial <- data.frame(subject = "s1", trial = 1, bin = 1)
  expect_error(
    ptl(trial, column = "missing", method = "shuffle"),
    "column not found"
  )
  expect_error(
    ptl(trial, column = "bin", method = "circular_shift", order_by = "missing"),
    "order_by column"
  )
})

test_that("residualize_trials skips groups with failed fits", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 4),
    x = rep(c(1, 1, 2, 2), 2),
    y = c(0, 1, 0, 1, 0, 1, 0, 1)
  )
  out <- residualize_trials(dat, response = "y", terms = "x", by = "subject", family = "binomial")
  expect_true(any(is.finite(out$y_resid)))
})

test_that("aggregate_by_bin handles empty bin grids", {
  dat <- data.frame(subject = character(0), bin = numeric(0), value = numeric(0))
  out <- aggregate_by_bin(dat, value = "value", bin = "bin", by = "subject", bins = numeric(0), complete = TRUE)
  expect_equal(nrow(out), 0L)
})

test_that("group_tfr_test validates replicate and cluster arguments", {
  x <- make_group_sine()
  expect_error(
    group_tfr_test(x, fs = 24, freqs = 6, nrep = 0),
    "positive scalar"
  )
  expect_error(
    group_tfr_test(x, fs = 24, freqs = 6, cluster_threshold = 0),
    "cluster_threshold"
  )
})

test_that("group_tfr resolves default and per-frequency bandwidth", {
  x <- make_group_sine(n_bin = 90)
  auto <- group_tfr(x, fs = 24, freqs = c(4, 6, 8), bandwidth = NULL, edge = 2)
  expect_equal(dim(auto$observed$map), c(3L, 90 - 4))

  custom <- group_tfr(x, fs = 24, freqs = c(4, 6), bandwidth = c(0.5, 0.75), edge = 2)
  expect_equal(dim(custom$observed$map), c(2L, 90 - 4))

  expect_error(group_tfr(x, fs = 24, freqs = 6, bandwidth = c(0.5, 0.75)), "same length as freqs")
  expect_error(group_tfr(x, fs = 24, freqs = c(0, 6)), "Nyquist")
})

test_that("group_tfr_test supports p_adjust=none without clustering", {
  x <- make_group_sine(n_subject = 4, n_bin = 72)
  res <- group_tfr_test(
    x,
    fs = 24,
    freqs = c(4, 6, 8),
    bandwidth = 1,
    null = "shuffle_labels",
    nrep = 8,
    p_adjust = "none",
    cluster = FALSE,
    edge = 2,
    seed = 44
  )
  expect_equal(res$stats$p_adj, res$stats$p)
  expect_equal(res$stats$clusters, list())
})

test_that("group_phase_difference validates matching subject grids", {
  x <- make_group_sine(n_subject = 3, n_bin = 60)
  y <- x[, -1, drop = FALSE]
  expect_error(
    group_phase_difference(x, y, fs = 24, freq = 6),
    "matching bins"
  )
})
