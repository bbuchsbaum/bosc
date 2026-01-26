test_that("oscore_tidy summarizes a single oscore result", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  res <- oscillation_score(sig, fs = fs, flim = c(1, 20), warnings = FALSE)
  df <- oscore_tidy(res)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_true(all(c("oscore", "fosc", "fmin", "fmax") %in% names(df)))
})

test_that("oscore_tidy binds a list of oscore results", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  res1 <- oscillation_score(sig, fs = fs, flim = c(1, 20), warnings = FALSE)
  res2 <- oscillation_score(sig, fs = fs, flim = c(5, 30), warnings = FALSE)
  df <- oscore_tidy(list(res1, res2))
  expect_equal(nrow(df), 2)
})

test_that("oscore_spectrum extracts spectrum from oscore output", {
  fs <- 200
  sig <- numeric(400)
  sig[seq(50, 350, 50)] <- 1
  res <- oscillation_score(sig, fs = fs, flim = c(1, 20), warnings = FALSE)
  sp <- oscore_spectrum(res)
  expect_s3_class(sp, "data.frame")
  expect_true(all(c("freq", "power") %in% names(sp)))
  expect_equal(nrow(sp), length(res$freqs))
})

test_that("clusters_tidy and clusters_pixels return consistent sizes", {
  mat <- matrix(rnorm(100), nrow = 10)
  mat[3:4, 3:4] <- mat[3:4, 3:4] + 3
  clusters <- detect_clusters(mat, threshold = 2, method = "extract")
  df_sum <- clusters_tidy(clusters)
  df_px <- clusters_pixels(clusters)
  expect_equal(nrow(df_sum), length(clusters))
  expect_equal(nrow(df_px), sum(vapply(clusters, function(cl) length(cl$Zscores), integer(1))))
})

# --- oscore_tidy additional branches ---

test_that("oscore_tidy errors on NULL input", {
  expect_error(oscore_tidy(NULL), "x cannot be NULL")
})

test_that("oscore_tidy errors on invalid input", {
  expect_error(oscore_tidy("bad"), "x must be an oscillation score result")
})

test_that("oscore_tidy includes z/pval/significant when present", {
  res <- list(oscore = 1.5, fosc = 10, flim = c(1, 20),
              z = 2.0, pval = 0.02, significant = TRUE)
  df <- oscore_tidy(res)
  expect_true(all(c("z", "pval", "significant") %in% names(df)))
  expect_equal(df$z, 2.0)
})

test_that("oscore_tidy handles result without flim", {
  res <- list(oscore = 1.5, fosc = 10)
  df <- oscore_tidy(res)
  expect_true(is.na(df$fmin))
  expect_true(is.na(df$fmax))
})

# --- oscore_spectrum additional branches ---

test_that("oscore_spectrum errors on NULL input", {
  expect_error(oscore_spectrum(NULL), "x cannot be NULL")
})

test_that("oscore_spectrum extracts fxx/spectrum fields", {
  res <- list(fxx = c(1, 2, 3), spectrum = c(0.1, 0.5, 0.2))
  df <- oscore_spectrum(res)
  expect_equal(nrow(df), 3)
  expect_true(all(c("freq", "power") %in% names(df)))
})

test_that("oscore_spectrum handles list of results with id column", {
  res1 <- list(fxx = c(1, 2), spectrum = c(0.1, 0.2))
  res2 <- list(fxx = c(3, 4), spectrum = c(0.3, 0.4))
  df <- oscore_spectrum(list(res1, res2))
  expect_true("id" %in% names(df))
  expect_equal(nrow(df), 4)
  expect_equal(sort(unique(df$id)), c(1, 2))
})

# --- clusters_tidy additional branches ---

test_that("clusters_tidy returns empty data.frame for NULL input", {
  df <- clusters_tidy(NULL)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0)
})

test_that("clusters_tidy returns empty data.frame for empty list", {
  df <- clusters_tidy(list())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0)
})

test_that("clusters_tidy handles cluster with NULL pAdj", {
  cl <- list(list(
    blobID = 1, times = c(1, 2), freqs = c(1, 2),
    Zscores = c(3, 4), CoMtime = 1.5, CoMfreq = 1.5,
    peakZscore = 4, sumZscore = 7, pAdj = NULL, peakpAdj = NULL
  ))
  df <- clusters_tidy(cl)
  expect_equal(nrow(df), 1)
  expect_true(is.na(df$minpAdj))
})

# --- clusters_pixels additional branches ---

test_that("clusters_pixels returns empty data.frame for NULL input", {
  df <- clusters_pixels(NULL)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0)
})

