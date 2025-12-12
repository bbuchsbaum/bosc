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
