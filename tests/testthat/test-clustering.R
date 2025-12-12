test_that("extract_clusters labels simple blobs", {
  mat <- matrix(0, nrow = 5, ncol = 5)
  mat[2, 2] <- 1
  mat[4, 4] <- 1
  lab <- extract_clusters(mat, thres = 0.5)
  labs <- setdiff(unique(lab), 0)
  expect_length(labs, 2)
})

test_that("detect_clusters finds a cluster and sumZ is defined", {
  mat <- matrix(0, nrow = 5, ncol = 5)
  mat[3, 3] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, split_clusters = FALSE)
  expect_length(res, 1)
  expect_true(!is.null(res[[1]]$sumZscore))
})

test_that("detect_clusters finds two separated clusters", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[2, 2] <- 5
  mat[8, 8] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, split_clusters = FALSE)
  expect_length(res, 2)
})
