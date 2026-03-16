test_that("align_time_bins maps values to nearest target bins", {
  dat <- data.frame(
    time = c(199.8, 233.7, 267.0, NA_real_)
  )

  out <- align_time_bins(
    dat,
    time = "time",
    bins = c(200, 233.33, 266.67),
    out = "time_aligned"
  )

  expect_equal(out$time_aligned[1], 200)
  expect_equal(out$time_aligned[2], 233.33)
  expect_equal(out$time_aligned[3], 266.67)
  expect_true(is.na(out$time_aligned[4]))
})

test_that("align_time_bins respects tolerance", {
  dat <- data.frame(time = c(200.2, 240))

  out <- align_time_bins(
    dat,
    time = "time",
    bins = c(200, 233.33),
    tolerance = 2
  )

  expect_equal(out$time_bin[1], 200)
  expect_true(is.na(out$time_bin[2]))
})

test_that("residualize_trials removes nuisance terms within groups", {
  x <- c(-2, -1, 0, 1, 2)
  signal <- c(2, -1, -2, -1, 2)
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 5),
    x = rep(x, 2),
    signal = rep(signal, 2)
  )
  dat$y <- c(10 + 3 * x + signal, -4 + 3 * x + signal)

  out <- residualize_trials(
    dat,
    response = "y",
    terms = "x",
    by = "subject",
    out = "y_resid"
  )

  expect_equal(out$y_resid, dat$signal, tolerance = 1e-10)
})

test_that("residualize_trials copies response when no nuisance terms are given", {
  dat <- data.frame(y = c(1, 2, NA_real_))
  out <- residualize_trials(dat, response = "y", terms = NULL)
  expect_equal(out$y_resid, dat$y)
})

test_that("aggregate_by_bin returns complete group-by-bin grid with counts", {
  dat <- data.frame(
    subject = c("s1", "s1", "s1", "s2"),
    bin = c(1, 1, 3, 2),
    value = c(10, 14, 7, 5)
  )

  out <- aggregate_by_bin(
    dat,
    value = "value",
    bin = "bin",
    by = "subject",
    bins = 1:3,
    complete = TRUE
  )

  expect_equal(nrow(out), 6L)
  s1b1 <- out[out$subject == "s1" & out$bin == 1, ]
  s1b2 <- out[out$subject == "s1" & out$bin == 2, ]
  s2b2 <- out[out$subject == "s2" & out$bin == 2, ]

  expect_equal(s1b1$value, 12)
  expect_equal(s1b1$n, 2)
  expect_equal(s1b1$n_nonmissing, 2)
  expect_true(is.na(s1b2$value))
  expect_equal(s1b2$n, 0)
  expect_equal(s2b2$value, 5)
})

test_that("aggregate_by_bin can retain only observed bins", {
  dat <- data.frame(
    subject = c("s1", "s1", "s2"),
    bin = c(1, 3, 2),
    value = c(10, 7, 5)
  )

  out <- aggregate_by_bin(
    dat,
    value = "value",
    bin = "bin",
    by = "subject",
    complete = FALSE
  )

  expect_equal(nrow(out), 3L)
  expect_false(any(out$subject == "s1" & out$bin == 2))
})

test_that("grouped-data helpers validate required columns", {
  dat <- data.frame(x = 1:3)

  expect_error(align_time_bins(dat, time = "missing", bins = c(1, 2)), "time column")
  expect_error(residualize_trials(dat, response = "missing", terms = "x"), "response column")
  expect_error(aggregate_by_bin(dat, value = "missing", bin = "x"), "value column")
})
