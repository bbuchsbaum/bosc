manual_group_sine <- function(n_subject = 5, n_bin = 96, fs = 24, freq = 6) {
  t <- seq(0, by = 1 / fs, length.out = n_bin)
  x <- vapply(seq_len(n_subject), function(i) {
    sin(2 * pi * freq * t + (i - 1) * pi / 7)
  }, numeric(n_bin))
  t(x)
}

test_that("residualize_trials matches manual grouped linear-model residuals", {
  x <- rep(c(-2, -1, 0, 1, 2), 2)
  z <- rep(c(1, -1, 1, -1, 0), 2)
  subject <- rep(c("s1", "s2"), each = 5)
  y <- c(5 + 2 * x[1:5] - 3 * z[1:5] + c(1, 0, -1, 0, 1),
         -2 + 2 * x[6:10] - 3 * z[6:10] + c(-1, 1, 0, 1, -1))
  dat <- data.frame(subject = subject, x = x, z = z, y = y)

  out <- residualize_trials(dat, response = "y", terms = c("x", "z"), by = "subject")

  manual <- unsplit(lapply(split(dat, dat$subject), function(df) {
    stats::residuals(stats::lm(y ~ x + z, data = df))
  }), dat$subject)

  expect_equal(unname(out$y_resid), unname(manual), tolerance = 1e-12)

  by(out, out$subject, function(df) {
    expect_true(abs(sum(df$y_resid)) < 1e-10)
    expect_true(abs(sum(df$y_resid * df$x)) < 1e-10)
    expect_true(abs(sum(df$y_resid * df$z)) < 1e-10)
  })
})

test_that("residualize_trials matches manual grouped binomial glm residuals", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 8),
    task = rep(c(-0.5, 0.5), times = 8),
    trial = rep(1:8, times = 2)
  )
  dat$y <- c(0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1)

  out <- residualize_trials(
    dat,
    response = "y",
    terms = c("task", "trial"),
    by = "subject",
    family = "binomial",
    type = "response"
  )

  manual <- unsplit(lapply(split(dat, dat$subject), function(df) {
    stats::residuals(
      stats::glm(y ~ task + trial, data = df, family = stats::binomial()),
      type = "response"
    )
  }), dat$subject)

  expect_equal(unname(out$y_resid), unname(manual), tolerance = 1e-12)
})

test_that("aggregate_by_bin is invariant to input row order", {
  dat <- data.frame(
    subject = c("s1", "s1", "s2", "s2", "s1", "s2"),
    cond = c("a", "a", "a", "a", "b", "b"),
    bin = c(1, 2, 1, 2, 1, 1),
    value = c(10, 20, 11, 21, 30, 31)
  )

  ref <- aggregate_by_bin(dat, value = "value", bin = "bin", by = c("subject", "cond"), bins = 1:2)
  perm <- aggregate_by_bin(dat[c(6, 3, 1, 5, 2, 4), ], value = "value", bin = "bin", by = c("subject", "cond"), bins = 1:2)

  expect_equal(ref, perm)
})

test_that("make_grouped_series reshapes trial data and drops incomplete groups", {
  dat <- data.frame(
    subject = c("s1", "s1", "s1", "s2", "s2"),
    bin = c(1, 2, 3, 1, 3),
    value = c(10, 20, 30, 11, 31)
  )

  out <- make_grouped_series(
    dat,
    value = "value",
    bin = "bin",
    subject = "subject",
    bins = 1:3,
    incomplete = "drop"
  )

  expect_equal(out$subjects, "s1")
  expect_equal(as.numeric(out$data[1, ]), c(10, 20, 30))
  expect_equal(out$bins, 1:3)
})

test_that("permute_trial_labels preserves within-group label multisets", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 5),
    trial = rep(1:5, times = 2),
    bin = c(1, 2, 3, 4, 5, 5, 4, 3, 2, 1)
  )

  shuf <- permute_trial_labels(dat, column = "bin", by = "subject", method = "shuffle", seed = 1)
  by(seq_len(nrow(dat)), dat$subject, function(idx) {
    expect_equal(sort(shuf$bin[idx]), sort(dat$bin[idx]))
  })

  shift <- permute_trial_labels(
    dat,
    column = "bin",
    by = "subject",
    method = "circular_shift",
    order_by = "trial",
    seed = 3
  )
  by(seq_len(nrow(dat)), dat$subject, function(idx) {
    expect_equal(sort(shift$bin[idx]), sort(dat$bin[idx]))
  })
})

test_that("group_spectrum matches manual per-subject spectral computation", {
  x <- manual_group_sine()
  preprocess <- getFromNamespace("preprocess_group_series", "bosc")

  res <- group_spectrum(x, fs = 24, flim = c(3, 10), detrend = TRUE, taper = "hann")

  manual_specs <- lapply(seq_len(nrow(x)), function(i) {
    spectral_peak(preprocess(x[i, ], detrend = TRUE), fs = 24, flim = c(3, 10), taper = "hann")
  })
  manual_power <- do.call(rbind, lapply(manual_specs, `[[`, "spectrum"))
  manual_freq <- manual_specs[[1]]$fxx
  manual_subject_peak <- vapply(manual_specs, `[[`, numeric(1), "freq")

  expect_equal(res$observed$freq, manual_freq, tolerance = 1e-12)
  expect_equal(res$observed$per_subject, manual_power, tolerance = 1e-12)
  expect_equal(res$observed$power, colMeans(manual_power), tolerance = 1e-12)
  expect_equal(res$observed$subject_peak, manual_subject_peak, tolerance = 1e-12)
})

test_that("group_spectrum is invariant to subject order and removable linear trends", {
  x <- manual_group_sine(n_subject = 6, n_bin = 100, fs = 25, freq = 5)
  t <- seq_len(ncol(x))
  trend <- outer(seq_len(nrow(x)), t, function(i, tt) 3 * i - 0.2 * tt)
  x_trend <- x + trend

  ref <- group_spectrum(x, fs = 25, flim = c(3, 8), detrend = TRUE, taper = "hann")
  perm <- group_spectrum(x[c(6, 2, 4, 1, 5, 3), ], fs = 25, flim = c(3, 8), detrend = TRUE, taper = "hann")
  det <- group_spectrum(x_trend, fs = 25, flim = c(3, 8), detrend = TRUE, taper = "hann")

  expect_equal(ref$observed$power, perm$observed$power, tolerance = 1e-12)
  expect_equal(ref$observed$power, det$observed$power, tolerance = 1e-10)
  expect_equal(ref$observed$peak, det$observed$peak, tolerance = 1e-10)
})

test_that("group_spectrum supports higher-order polynomial detrending", {
  x <- manual_group_sine(n_subject = 4, n_bin = 96, fs = 24, freq = 6)
  t <- seq(-1, 1, length.out = ncol(x))
  quad <- outer(seq_len(nrow(x)), t^2, function(i, tt) 4 * i * tt)
  x_quad <- x + quad

  ref <- group_spectrum(x, fs = 24, flim = c(3, 10), detrend = TRUE, detrend_order = 2, taper = "hann")
  det <- group_spectrum(x_quad, fs = 24, flim = c(3, 10), detrend = TRUE, detrend_order = 2, taper = "hann")

  expect_equal(ref$observed$power, det$observed$power, tolerance = 1e-8)
  expect_equal(ref$observed$peak, det$observed$peak, tolerance = 1e-8)
})

test_that("group_spectrum pad_to yields exact FFT-bin frequencies", {
  x <- manual_group_sine(n_subject = 3, n_bin = 28, fs = 30, freq = 7)

  res <- group_spectrum(
    x,
    fs = 30,
    flim = c(3, 14),
    detrend = TRUE,
    detrend_order = 2,
    pad_to = 30,
    taper = "hann"
  )

  expect_equal(res$observed$freq, 3:14, tolerance = 1e-12)
  expect_equal(res$observed$peak, 7, tolerance = 1e-12)
})

test_that("group_spectrum raw_power matches direct FFT power bins", {
  x <- manual_group_sine(n_subject = 2, n_bin = 28, fs = 30, freq = 7)

  res <- group_spectrum(
    x,
    fs = 30,
    flim = c(3, 14),
    detrend = FALSE,
    pad_to = 30,
    spectrum = "raw_power",
    taper = "hann"
  )

  manual <- lapply(seq_len(nrow(x)), function(i) {
    xi <- x[i, ] * 0.5 * (1 - cos(2 * pi * seq(0, 27) / 27))
    xi <- c(xi, 0, 0)
    fft_vals <- fft(xi)
    spec <- Mod(fft_vals)^2
    spec[4:15]
  })
  manual <- do.call(rbind, manual)

  expect_equal(res$observed$freq, 3:14, tolerance = 1e-12)
  expect_equal(res$observed$per_subject, manual, tolerance = 1e-10)
  expect_equal(res$observed$power, colMeans(manual), tolerance = 1e-10)
})

test_that("group_spectrum_test null draws are exactly reproducible from internals", {
  x <- manual_group_sine(n_subject = 4, n_bin = 72, fs = 24, freq = 6)
  gen_null <- getFromNamespace("generate_group_spectrum_null", "bosc")
  comp_spec <- getFromNamespace("compute_group_spectrum", "bosc")

  res <- group_spectrum_test(
    x,
    fs = 24,
    flim = c(3, 10),
    null = "circular_shift",
    nrep = 6,
    detrend = TRUE,
    taper = "hann",
    seed = 123
  )

  set.seed(123)
  manual_null <- array(NA_real_, dim = c(6, length(res$observed$freq)))
  for (i in 1:6) {
    draw <- gen_null(x, method = "circular_shift")
    manual_null[i, ] <- comp_spec(
      draw,
      fs = 24,
      flim = c(3, 10),
      detrend = TRUE,
      detrend_order = 1,
      taper = "hann",
      fcor = FALSE
    )$power
  }

  expect_equal(res$null$power, manual_null, tolerance = 1e-12)
})

test_that("group_spectrum_test_trials is reproducible from explicit trial shuffles", {
  x <- manual_group_sine(n_subject = 4, n_bin = 72, fs = 24, freq = 6)
  trial_dat <- data.frame(
    subject = rep(paste0("s", seq_len(nrow(x))), each = ncol(x)),
    bin = rep(seq_len(ncol(x)), times = nrow(x)),
    value = as.numeric(t(x))
  )

  trial_res <- group_spectrum_test_trials(
    trial_dat,
    value = "value",
    bin = "bin",
    subject = "subject",
    fs = 24,
    bins = seq_len(ncol(x)),
    flim = c(3, 10),
    null = "shuffle_bins",
    nrep = 5,
    detrend = TRUE,
    detrend_order = 1,
    taper = "hann",
    seed = 123
  )

  comp_spec <- getFromNamespace("compute_group_spectrum", "bosc")
  set.seed(123)
  manual_null <- array(NA_real_, dim = c(5, length(trial_res$observed$freq)))
  for (i in 1:5) {
    draw <- permute_trial_labels(trial_dat, column = "bin", by = "subject", method = "shuffle")
    series <- make_grouped_series(
      draw,
      value = "value",
      bin = "bin",
      subject = "subject",
      bins = seq_len(ncol(x))
    )
    manual_null[i, ] <- comp_spec(
      series$data,
      fs = 24,
      flim = c(3, 10),
      detrend = TRUE,
      detrend_order = 1,
      taper = "hann",
      fcor = FALSE
    )$power
  }

  expect_equal(trial_res$null$power, manual_null, tolerance = 1e-12)
})

test_that("group_tfr matches manual narrowband Hilbert amplitudes for one frequency", {
  x <- manual_group_sine(n_subject = 4, n_bin = 80, fs = 20, freq = 4)
  preprocess <- getFromNamespace("preprocess_group_series", "bosc")

  res <- group_tfr(x, fs = 20, freqs = 4, bandwidth = 0.75, detrend = TRUE, reflect = FALSE, edge = 0)

  manual <- vapply(seq_len(nrow(x)), function(i) {
    nb <- narrowband_hilbert(preprocess(x[i, ], detrend = TRUE), fs = 20, freqlim = c(3.25, 4.75))
    Mod(nb$analytic)
  }, numeric(ncol(x)))
  manual <- t(manual)

  expect_equal(dim(res$observed$per_subject), c(nrow(x), 1L, ncol(x)))
  expect_equal(as.numeric(res$observed$per_subject[, 1, ]), as.numeric(manual), tolerance = 1e-10)
  expect_equal(as.numeric(res$observed$map[1, ]), colMeans(manual), tolerance = 1e-10)
})

test_that("group_tfr is invariant to subject order", {
  x <- manual_group_sine(n_subject = 5, n_bin = 90, fs = 30, freq = 6)
  ref <- group_tfr(x, fs = 30, freqs = c(4, 6, 8), bandwidth = 1, reflect = TRUE, edge = 3)
  perm <- group_tfr(x[c(5, 1, 4, 2, 3), ], fs = 30, freqs = c(4, 6, 8), bandwidth = 1, reflect = TRUE, edge = 3)

  expect_equal(ref$observed$map, perm$observed$map, tolerance = 1e-12)
})

test_that("group_tfr_test null maps are reproducible from internals", {
  x <- manual_group_sine(n_subject = 3, n_bin = 60, fs = 20, freq = 4)
  gen_null <- getFromNamespace("generate_group_spectrum_null", "bosc")
  comp_tfr <- getFromNamespace("compute_group_tfr", "bosc")

  res <- group_tfr_test(
    x,
    fs = 20,
    freqs = c(3, 4),
    bandwidth = 0.5,
    null = "shuffle_labels",
    nrep = 4,
    reflect = FALSE,
    edge = 2,
    seed = 222,
    cluster = FALSE
  )

  set.seed(222)
  manual_null <- array(NA_real_, dim = c(4, 2, 56))
  for (i in 1:4) {
    draw <- gen_null(x, method = "shuffle_labels")
    manual_null[i, , ] <- comp_tfr(
      draw,
      fs = 20,
      freqs = c(3, 4),
      bandwidth = c(0.5, 0.5),
      detrend = TRUE,
      reflect = FALSE,
      edge = 2,
      filtorder = 2,
      demean = TRUE,
      tol = 1e2
    )$map
  }

  expect_equal(res$null$map, manual_null, tolerance = 1e-10)
})

test_that("group_phase_consistency returns perfect concentration for identical subjects", {
  t <- seq(0, by = 1 / 20, length.out = 80)
  sig <- sin(2 * pi * 4 * t)
  x <- matrix(rep(sig, times = 5), nrow = 5, byrow = TRUE)

  res <- group_phase_consistency(x, fs = 20, freq = 4, bandwidth = 0.75, reflect = FALSE, edge = 0)

  expect_true(all(abs(res$stats$r - 1) < 1e-10))
  expect_true(all(res$stats$p <= 0.05))
})

test_that("group_phase_difference is zero for identical inputs and wrapped as specified", {
  pair <- manual_group_sine(n_subject = 4, n_bin = 72, fs = 24, freq = 6)
  same <- group_phase_difference(pair, pair, fs = 24, freq = 6, bandwidth = 0.75, reflect = FALSE, edge = 0)

  expect_true(all(abs(same$observed$phase_difference) < 1e-10))
  expect_true(all(abs(same$observed$mean_phase_difference) < 1e-10))
  expect_true(all(abs(same$stats$r - 1) < 1e-10))

  t <- seq(0, by = 1 / 24, length.out = 72)
  x1 <- matrix(rep(sin(2 * pi * 6 * t), times = 4), nrow = 4, byrow = TRUE)
  x2 <- matrix(rep(sin(2 * pi * 6 * t + pi / 2), times = 4), nrow = 4, byrow = TRUE)
  wrapped <- group_phase_difference(x1, x2, fs = 24, freq = 6, bandwidth = 0.75, mode = "wrapped_2pi")

  expected <- 3 * pi / 2
  actual <- circ_mean(wrapped$observed$mean_phase_difference)
  delta <- ((actual - expected) + pi) %% (2 * pi) - pi
  expect_true(abs(delta) < 0.2)
})
