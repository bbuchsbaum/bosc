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

