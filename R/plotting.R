# Declare global variables used in ggplot2 aes() to avoid R CMD check NOTEs
utils::globalVariables(c("freq", "power", "phase", "time", "value", "x", "y"))

#' Plot spectrum
#'
#' Convenience wrapper to visualize spectra returned by \code{spectral_peak}.
#'
#' @param freqs frequency axis.
#' @param spectrum power values.
#' @param peak optional peak frequency to highlight.
#' @return ggplot object.
#' @export
plot_spectrum <- function(freqs, spectrum, peak = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting.")
  }
  df <- data.frame(freq = freqs, power = spectrum)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = freq, y = power)) +
    ggplot2::geom_line(color = "#1B75BB") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "Frequency (Hz)", y = "Power") +
    ggplot2::theme_minimal()
  if (!is.null(peak) && is.finite(peak)) {
    p <- p + ggplot2::geom_vline(xintercept = peak, linetype = "dashed", color = "#0F4C81")
  }
  p
}

#' Plot phase histogram
#'
#' @param phases numeric vector of phases (radians).
#' @param nbins number of bins.
#' @return ggplot object.
#' @export
plot_phase_hist <- function(phases, nbins = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting.")
  }
  breaks <- seq(-pi, pi, length.out = nbins + 1)
  df <- data.frame(phase = phases)
  ggplot2::ggplot(df, ggplot2::aes(x = phase)) +
    ggplot2::geom_histogram(breaks = breaks, fill = "#1B75BB", color = "white", boundary = -pi) +
    ggplot2::scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                                labels = c("-pi", "-pi/2", "0", "pi/2", "pi")) +
    ggplot2::labs(x = "Phase (rad)", y = "Count") +
    ggplot2::theme_minimal()
}

#' Simple Box Plot by Group
#'
#' Creates narrow-form box plots showing quantile ranges (5th-95th percentile as
#' thin line, 25th-75th as thick line, median as point). Port of MATLAB
#' \code{simplebox}. For groups with fewer than 4 observations, individual points
#' are shown instead.
#'
#' @param labels numeric or factor vector of group labels.
#' @param data numeric vector of data values (same length as \code{labels}).
#' @param colors optional matrix of RGB colors (one row per unique label), or a
#'   single color vector applied to all groups. Values should be 0-1 scale.
#'   If NULL, uses a default blue palette.
#' @return ggplot object. Also returns (invisibly) a vector of median values per group.
#' @export
#' @examples
#' set.seed(1)
#' labels <- rep(1:3, each = 20)
#' data <- rnorm(60) + labels * 0.5
#' plot_simplebox(labels, data)
plot_simplebox <- function(labels, data, colors = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting.")
  }
  if (length(labels) != length(data)) {
    stop("Data and labels do not match")
  }

  df <- data.frame(label = labels, value = data)
  ulabels <- unique(labels)
  n_groups <- length(ulabels)

  # Set up colors
 if (is.null(colors)) {
    if (requireNamespace("RColorBrewer", quietly = TRUE) && n_groups <= 9) {
      pal <- RColorBrewer::brewer.pal(max(3, n_groups), "Set1")
      colors <- grDevices::col2rgb(pal[seq_len(n_groups)]) / 255
      colors <- t(colors)
    } else {
      colors <- matrix(rep(c(0.11, 0.46, 0.73), n_groups), ncol = 3, byrow = TRUE)
    }
  } else if (is.vector(colors) && length(colors) == 3) {
    colors <- matrix(rep(colors, n_groups), ncol = 3, byrow = TRUE)
  }

  # Compute quantiles per group
  stats_list <- lapply(seq_along(ulabels), function(i) {
    lb <- ulabels[i]
    datalb <- df$value[df$label == lb]
    col_hex <- grDevices::rgb(colors[i, 1], colors[i, 2], colors[i, 3])

    if (length(datalb) > 3) {
      q <- stats::quantile(datalb, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
      list(
        label = lb,
        q05 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q95 = q[5],
        color = col_hex,
        show_box = TRUE,
        points = NULL
      )
    } else {
      list(
        label = lb,
        q50 = if (length(datalb) > 0) stats::median(datalb) else NA,
        color = col_hex,
        show_box = FALSE,
        points = datalb
      )
    }
  })

  # Build plot
  p <- ggplot2::ggplot() + ggplot2::theme_minimal()

  for (s in stats_list) {
    if (s$show_box) {
      # 5th-95th percentile line (thin)
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = s$label, xend = s$label, y = s$q05, yend = s$q95),
        linewidth = 0.5, color = s$color
      )
      # 25th-75th percentile line (thick)
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = s$label, xend = s$label, y = s$q25, yend = s$q75),
        linewidth = 2, color = s$color
      )
      # Median point
      p <- p + ggplot2::geom_point(
        ggplot2::aes(x = s$label, y = s$q50),
        size = 3, color = s$color
      )
    } else if (!is.null(s$points) && length(s$points) > 0) {
      # Individual points for small groups
      pts_df <- data.frame(x = rep(s$label, length(s$points)), y = s$points)
      p <- p + ggplot2::geom_point(
        data = pts_df,
        ggplot2::aes(x = x, y = y),
        size = 2, color = s$color
      )
    }
  }

  p <- p + ggplot2::labs(x = "Group", y = "Value")

  # Return medians invisibly
  medians <- vapply(stats_list, function(s) {
    if (s$show_box) s$q50 else if (!is.null(s$points)) stats::median(s$points) else NA_real_
  }, numeric(1))
  names(medians) <- as.character(ulabels)

  attr(p, "medians") <- medians
  p
}

#' Plot Cluster Heatmap
#'
#' Visualize detected clusters as a heatmap with optional cluster boundary overlays.
#'
#' @param data matrix of values (e.g., Z-scores).
#' @param clusters optional list of cluster objects from \code{detect_clusters}.
#' @param freq_axis optional frequency axis for y-axis labels.
#' @param time_axis optional time axis for x-axis labels.
#' @param zlim optional length-2 vector for color scale limits.
#' @param palette color palette name (default "RdBu").
#' @return ggplot object.
#' @export
#' @examples
#' mat <- matrix(rnorm(400), nrow = 20)
#' mat[8:12, 8:12] <- mat[8:12, 8:12] + 3
#' plot_cluster_heatmap(mat)
plot_cluster_heatmap <- function(data, clusters = NULL,
                                  freq_axis = NULL, time_axis = NULL,
                                  zlim = NULL, palette = "RdBu") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting.")
  }

  if (is.null(freq_axis)) freq_axis <- seq_len(nrow(data))
  if (is.null(time_axis)) time_axis <- seq_len(ncol(data))

  # Create long-form data frame
  df <- expand.grid(time = time_axis, freq = freq_axis)
  df$value <- as.vector(t(data))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = freq, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed(ratio = diff(range(time_axis)) / diff(range(freq_axis)) * 0.5)

  # Color scale
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    pal_colors <- rev(RColorBrewer::brewer.pal(11, palette))
    if (!is.null(zlim)) {
      p <- p + ggplot2::scale_fill_gradientn(colors = pal_colors, limits = zlim,
                                              oob = scales::squish)
    } else {
      p <- p + ggplot2::scale_fill_gradientn(colors = pal_colors)
    }
  } else {
    p <- p + ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                            midpoint = 0)
  }

  # Add cluster boundaries if provided
  if (!is.null(clusters) && length(clusters) > 0) {
    for (cl in clusters) {
      if (!is.null(cl$times) && !is.null(cl$freqs)) {
        hull_df <- data.frame(time = cl$times, freq = cl$freqs)
        # Simple convex hull outline
        if (nrow(hull_df) >= 3) {
          hull_idx <- grDevices::chull(hull_df$time, hull_df$freq)
          hull_pts <- hull_df[c(hull_idx, hull_idx[1]), ]
          p <- p + ggplot2::geom_path(data = hull_pts,
                                       ggplot2::aes(x = time, y = freq),
                                       color = "black", linewidth = 0.8,
                                       inherit.aes = FALSE)
        }
      }
    }
  }

  p <- p +
    ggplot2::labs(x = "Time", y = "Frequency", fill = "Z") +
    ggplot2::theme_minimal()

  p
}
