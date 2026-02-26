test_that("extract_clusters labels simple blobs", {
  mat <- matrix(0, nrow = 5, ncol = 5)
  mat[2, 2] <- 1
  mat[4, 4] <- 1
  lab <- extract_clusters(mat, thres = 0.5)
  labs <- setdiff(unique(lab), 0)
  expect_length(labs, 2)
})

test_that("extract_clusters validates matrix input and tiny dimensions", {
  expect_error(extract_clusters(1:5, thres = NULL), "matrix")

  tiny <- matrix(1, nrow = 2, ncol = 2)
  lab <- extract_clusters(tiny, thres = NULL)
  expect_true(all(lab == 0L))
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

test_that("extract_clusters matches reference 8-connectivity on random grids", {
  ref_cc <- function(mat) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    dat <- mat
    dat[1, ] <- 0
    dat[nr, ] <- 0
    dat[, 1] <- 0
    dat[, nc] <- 0
    out <- matrix(0L, nrow = nr, ncol = nc)
    neigh <- as.matrix(expand.grid(dr = -1:1, dc = -1:1))
    neigh <- neigh[!(neigh[, 1] == 0 & neigh[, 2] == 0), , drop = FALSE]
    lbl <- 1L
    for (r in 2:(nr - 1)) {
      for (c in 2:(nc - 1)) {
        if (dat[r, c] != 0 && out[r, c] == 0L) {
          qr <- integer(nr * nc)
          qc <- integer(nr * nc)
          head <- 1L
          tail <- 1L
          qr[tail] <- r
          qc[tail] <- c
          out[r, c] <- lbl
          while (head <= tail) {
            rr <- qr[head]
            cc <- qc[head]
            head <- head + 1L
            for (k in seq_len(nrow(neigh))) {
              r2 <- rr + neigh[k, 1]
              c2 <- cc + neigh[k, 2]
              if (r2 >= 2 && r2 <= (nr - 1) && c2 >= 2 && c2 <= (nc - 1) &&
                  dat[r2, c2] != 0 && out[r2, c2] == 0L) {
                tail <- tail + 1L
                qr[tail] <- r2
                qc[tail] <- c2
                out[r2, c2] <- lbl
              }
            }
          }
          lbl <- lbl + 1L
        }
      }
    }
    out
  }

  set.seed(123)
  for (i in seq_len(40)) {
    p <- stats::runif(1, 0.05, 0.35)
    mat <- matrix(stats::rbinom(12 * 12, 1, p), nrow = 12)
    lab <- extract_clusters(mat, thres = NULL)
    ref <- ref_cc(mat)
    active <- which(mat != 0 & row(mat) > 1 & row(mat) < 12 & col(mat) > 1 & col(mat) < 12, arr.ind = TRUE)
    if (nrow(active) < 2) next
    same_lab <- outer(seq_len(nrow(active)), seq_len(nrow(active)), Vectorize(function(a, b) {
      lab[active[a, 1], active[a, 2]] == lab[active[b, 1], active[b, 2]]
    }))
    same_ref <- outer(seq_len(nrow(active)), seq_len(nrow(active)), Vectorize(function(a, b) {
      ref[active[a, 1], active[a, 2]] == ref[active[b, 1], active[b, 2]]
    }))
    expect_true(all(same_lab == same_ref))
  }
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

test_that("detect_clusters can filter clusters using qval", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[3, 3] <- 5
  mat[8, 8] <- 1.1

  all_res <- detect_clusters(
    mat, threshold = 1, method = "extract", smooth = NULL,
    qval = 0.05, filter_by_q = FALSE
  )
  filt_res <- detect_clusters(
    mat, threshold = 1, method = "extract", smooth = NULL,
    qval = 0.05, filter_by_q = TRUE
  )

  expect_length(all_res, 2)
  expect_length(filt_res, 1)
  expect_true(abs(filt_res[[1]]$peakZscore) >= 5)
})

test_that("detect_clusters validates q filtering args", {
  mat <- matrix(0, nrow = 5, ncol = 5)
  mat[3, 3] <- 2
  expect_error(
    detect_clusters(mat, threshold = 1, method = "extract", filter_by_q = TRUE),
    "requires qval"
  )
  expect_error(
    detect_clusters(mat, threshold = 1, method = "extract", qval = 2),
    "qval must be in"
  )
})

test_that("detect_clusters with split_clusters=TRUE attempts splitting", {
  mat <- matrix(0, nrow = 20, ncol = 20)
  # Create a wide cluster that might be split
  mat[5:7, 3:5] <- 5
  mat[5:7, 15:17] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL, split_clusters = TRUE)
  expect_true(length(res) >= 1)
})

test_that("detect_clusters split reassignment works for bridged peaks", {
  mat <- matrix(0, nrow = 30, ncol = 30)
  # Two strong lobes connected by a weak bridge so splitting creates missing pixels.
  mat[10:14, 6:10] <- 5
  mat[10:14, 20:24] <- 5
  mat[11:13, 11:19] <- 1.1

  res <- detect_clusters(
    mat,
    threshold = 1,
    method = "extract",
    split_clusters = TRUE,
    peak_height = seq(0.2, 0.8, by = 0.1),
    min_new_cluster_size = 0.05,
    smooth = NULL
  )
  expect_true(length(res) >= 2)
  expect_true(all(vapply(res, function(cl) length(cl$times) > 0, logical(1))))
})

test_that("detect_clusters supports smoothing and watershed paths", {
  skip_if_not_installed("imager")
  mat <- matrix(0, nrow = 20, ncol = 20)
  mat[8:10, 8:10] <- 5
  mat[12:14, 12:14] <- -4

  smooth_res <- detect_clusters(
    mat,
    threshold = 1,
    method = "extract",
    smooth = 1
  )
  expect_true(length(smooth_res) >= 1)

  ws_res <- detect_clusters(
    mat,
    threshold = 1,
    method = "watershed",
    smooth = NULL
  )
  expect_true(is.list(ws_res))
})

test_that("detect_clusters errors when imager is unavailable for smooth/watershed", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 4
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "imager")) return(FALSE)
      base::requireNamespace(package, ...)
    },
    .package = "base"
  )

  expect_error(
    detect_clusters(mat, threshold = 1, method = "extract", smooth = 1),
    "imager"
  )
  expect_error(
    detect_clusters(mat, threshold = 1, method = "watershed", smooth = NULL),
    "imager"
  )
})

test_that("detect_clusters handles single-pixel cluster", {
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 5
  res <- detect_clusters(mat, threshold = 1, method = "extract", smooth = NULL)
  expect_length(res, 1)
  expect_length(res[[1]]$Zscores, 1)
})
