test_that("plot_spectrum creates a ggplot", {
  skip_if_not_installed("ggplot2")
  freqs <- seq(1, 10, by = 1)
  spec <- freqs * 0 + 1
  p <- plot_spectrum(freqs, spec, peak = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_phase_hist creates a ggplot", {
  skip_if_not_installed("ggplot2")
  phases <- seq(-pi, pi, length.out = 20)
  p <- plot_phase_hist(phases, nbins = 6)
  expect_s3_class(p, "ggplot")
})

# --- plot_spectrum additional branches ---

test_that("plot_spectrum without peak omits vline", {
  skip_if_not_installed("ggplot2")
  freqs <- seq(1, 10, by = 1)
  spec <- rep(1, length(freqs))
  p <- plot_spectrum(freqs, spec, peak = NULL)
  expect_s3_class(p, "ggplot")
  # Should have no geom_vline layer (only geom_line)
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomVline" %in% layer_classes)
})

test_that("plot_spectrum with NaN peak omits vline", {
  skip_if_not_installed("ggplot2")
  freqs <- seq(1, 10, by = 1)
  spec <- rep(1, length(freqs))
  p <- plot_spectrum(freqs, spec, peak = NaN)
  expect_s3_class(p, "ggplot")
  layer_classes <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomVline" %in% layer_classes)
})

# --- plot_phase_hist default nbins ---

test_that("plot_phase_hist uses default nbins = 12", {
  skip_if_not_installed("ggplot2")
  phases <- seq(-pi, pi, length.out = 50)
  p <- plot_phase_hist(phases)
  expect_s3_class(p, "ggplot")
})

# --- plot_simplebox ---

test_that("plot_simplebox creates box plot for large groups", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  labels <- rep(1:3, each = 20)
  data <- rnorm(60) + labels * 0.5
  p <- plot_simplebox(labels, data)
  expect_s3_class(p, "ggplot")
  medians <- attr(p, "medians")
  expect_length(medians, 3)
  expect_true(all(is.finite(medians)))
})

test_that("plot_simplebox shows points for small groups", {
  skip_if_not_installed("ggplot2")
  labels <- c(1, 1, 2, 2)
  data <- c(1, 2, 3, 4)
  p <- plot_simplebox(labels, data)
  expect_s3_class(p, "ggplot")
})

test_that("plot_simplebox accepts single-color vector", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  labels <- rep(1:2, each = 10)
  data <- rnorm(20)
  p <- plot_simplebox(labels, data, colors = c(0.5, 0.2, 0.8))
  expect_s3_class(p, "ggplot")
})

test_that("plot_simplebox errors on mismatched labels/data", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_simplebox(1:3, 1:5), "Data and labels do not match")
})

# --- plot_cluster_heatmap ---

test_that("plot_cluster_heatmap creates basic heatmap", {
  skip_if_not_installed("ggplot2")
  mat <- matrix(rnorm(100), nrow = 10)
  p <- plot_cluster_heatmap(mat)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cluster_heatmap draws cluster hulls", {
  skip_if_not_installed("ggplot2")
  mat <- matrix(rnorm(100), nrow = 10)
  mat[3:5, 3:5] <- mat[3:5, 3:5] + 5
  clusters <- detect_clusters(mat, threshold = 2, method = "extract", smooth = NULL)
  p <- plot_cluster_heatmap(mat, clusters = clusters)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cluster_heatmap respects zlim", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")
  mat <- matrix(rnorm(100), nrow = 10)
  p <- plot_cluster_heatmap(mat, zlim = c(-2, 2))
  expect_s3_class(p, "ggplot")
})

test_that("plot_cluster_heatmap skips hull for < 3 points", {
  skip_if_not_installed("ggplot2")
  mat <- matrix(0, nrow = 10, ncol = 10)
  mat[5, 5] <- 5
  # single-pixel cluster: hull should be skipped gracefully
  clusters <- detect_clusters(mat, threshold = 2, method = "extract", smooth = NULL)
  p <- plot_cluster_heatmap(mat, clusters = clusters)
  expect_s3_class(p, "ggplot")
})
