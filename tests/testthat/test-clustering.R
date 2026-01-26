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

# --- extract_clusters with NULL threshold (binary input) ---

test_that("extract_clusters with thres=NULL uses binary input directly", {
  mat <- matrix(0L, nrow = 5, ncol = 5)
  mat[3, 3] <- 1L
  lab <- extract_clusters(mat, thres = NULL)
  labs <- setdiff(unique(lab), 0)
  expect_length(labs, 1)
})

# --- detect_clusters additional branches ---

test_that("detect_clusters returns empty list for all-zero matrix", {
  mat <- matrix(0, nrow = 5, ncol = 5)
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL)
  expect_length(res, 0)
})

test_that("detect_clusters finds negative clusters", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- -5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL)
  expect_true(length(res) >= 1)
  # At least one cluster should have a negative blobID
  blob_ids <- vapply(res, function(cl) cl$blobID, numeric(1))
  expect_true(any(blob_ids < 0))
})

test_that("detect_clusters uses data2 for peak extraction", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 5
  data2 <- matrix(100, nrow = 10, ncol = 10)
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, data2 = data2)
  expect_equal(res[[1]]$peakZscore, 100)
})

test_that("detect_clusters uses dataRef for reference stats", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 5
  ref <- matrix(42, nrow = 10, ncol = 10)
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, dataRef = ref)
  expect_true(!is.null(res[[1]]$peakZscoreRef))
  expect_equal(res[[1]]$peakZscoreRef, 42)
})

test_that("detect_clusters computes pAdj with qval", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, qval = 0.05)
  expect_true(!is.na(res[[1]]$peakpAdj))
})

test_that("detect_clusters with split_clusters=TRUE attempts splitting", {
  mat <- matrix(0, nrow = 20, ncol = 20)
  # Create a wide cluster that might be split
  mat[5:7, 3:5] <- 5
  mat[5:7, 15:17] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, split_clusters = TRUE)
  expect_true(length(res) >= 1)
})

test_that("detect_clusters handles single-pixel cluster", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL)
  expect_length(res, 1)
  expect_length(res[[1]]$Zscores, 1)
})
