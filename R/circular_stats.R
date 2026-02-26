#' Circular Mean Direction
#'
#' Computes the mean direction for circular data. Port of MATLAB \code{circ_mean}
#' from the Circular Statistics Toolbox (Berens, 2009).
#'
#' @param alpha numeric vector or matrix of angles in radians.
#' @param w optional weights (same dimensions as \code{alpha}). If NULL, uniform weights.
#' @param dim integer; dimension along which to compute mean (default 1).
#' @return If \code{alpha} is a vector, returns a single mean direction. If matrix,
#'   returns a vector of means computed along \code{dim}.
#' @references
#'   Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics.
#'   Journal of Statistical Software, 31(10), 1-21.
#'
#'   Fisher NI (1993). Statistical Analysis of Circular Data. Cambridge University Press.
#' @seealso \code{\link{circ_r}}, \code{\link{circ_vtest}}
#' @export
#' @examples
#' angles <- c(0, pi/4, pi/2)
#' circ_mean(angles)
#'
#' # Weighted mean
#' circ_mean(angles, w = c(1, 2, 1))
circ_mean <- function(alpha, w = NULL, dim = 1) {
  alpha <- as.array(alpha)
  dims <- dim(alpha)


  if (is.null(w)) {
    w <- array(1, dim = dims)
  } else {
    w <- as.array(w)
    if (!all(dim(w) == dims)) {
      stop("Input dimensions do not match")
    }
  }

  # For vectors, treat as 1D

if (length(dims) == 0 || (length(dims) == 1 && dims[1] == length(alpha))) {
    r <- sum(w * exp(1i * alpha))
    return(Arg(r))
  }

  # For matrices/arrays, apply along dimension
  r <- apply(w * exp(1i * alpha), setdiff(seq_along(dims), dim), sum)
  Arg(r)
}

#' Mean Resultant Vector Length
#'
#' Computes the mean resultant vector length for circular data, a measure of
#' concentration (0 = uniform, 1 = perfectly aligned). Port of MATLAB \code{circ_r}
#' from the Circular Statistics Toolbox (Berens, 2009).
#'
#' @param alpha numeric vector or matrix of angles in radians.
#' @param w optional weights (same dimensions as \code{alpha}). If NULL, uniform weights.
#' @param d optional spacing of bin centers for binned data (radians). If provided,
#'   applies correction factor for bias in estimation of r.
#' @param dim integer; dimension along which to compute (default 1).
#' @return Mean resultant length (scalar for vector input, vector/array otherwise).
#' @references
#'   Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics.
#'   Journal of Statistical Software, 31(10), 1-21.
#'
#'   Zar JH (2010). Biostatistical Analysis. 5th ed. Prentice Hall. (Equation 26.16)
#' @seealso \code{\link{circ_mean}}, \code{\link{circ_vtest}}
#' @export
#' @examples
#' # Highly concentrated data
#' angles <- c(0, 0.1, -0.1, 0.05)
#' circ_r(angles)  # Should be close to 1
#'
#' # Uniformly distributed data
#' angles_uniform <- seq(0, 2*pi, length.out = 100)
#' circ_r(angles_uniform)  # Should be close to 0
circ_r <- function(alpha, w = NULL, d = 0, dim = 1) {
  alpha <- as.array(alpha)
  dims <- dim(alpha)

  if (is.null(w)) {
    w <- array(1, dim = dims)
  } else {
    w <- as.array(w)
    if (!all(dim(w) == dims)) {
      stop("Input dimensions do not match")
    }
  }

  # For vectors
  if (length(dims) == 0 || (length(dims) == 1 && dims[1] == length(alpha))) {
    r_complex <- sum(w * exp(1i * alpha))
    r <- abs(r_complex) / sum(w)
    # Apply binning correction if d != 0
    if (d != 0) {
      c_corr <- d / 2 / sin(d / 2)
      r <- c_corr * r
    }
    return(r)
  }

  # For matrices/arrays
  r_complex <- apply(w * exp(1i * alpha), setdiff(seq_along(dims), dim), sum)
  w_sum <- apply(w, setdiff(seq_along(dims), dim), sum)
  r <- abs(r_complex) / w_sum

  if (d != 0) {
    c_corr <- d / 2 / sin(d / 2)
    r <- c_corr * r
  }
  r
}

#' V-test for Non-uniformity with Specified Mean Direction
#'
#' Computes the V-test for circular data, testing whether the population is
#' non-uniformly distributed with a specified mean direction. Port of MATLAB
#' \code{circ_vtest} from the Circular Statistics Toolbox (Berens, 2009).
#'
#' The V-test is more powerful than the Rayleigh test when there is reason to
#' believe in a specific mean direction.
#'
#' @param alpha numeric vector of angles in radians.
#' @param dir suspected mean direction in radians.
#' @param w optional weights for binned data (same length as \code{alpha}).
#' @param d optional bin spacing for binned data (radians) for bias correction.
#' @return list with:
#'   \describe{
#'     \item{pval}{p-value of the V-test (one-tailed)}
#'     \item{v}{V statistic}
#'   }
#' @references
#'   Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics.
#'   Journal of Statistical Software, 31(10), 1-21.
#'
#'   Zar JH (2010). Biostatistical Analysis. 5th ed. Prentice Hall.
#' @seealso \code{\link{circ_mean}}, \code{\link{circ_r}}
#' @export
#' @examples
#' # Test if data clusters around 0
#' angles <- rnorm(50, mean = 0, sd = 0.5)
#' result <- circ_vtest(angles, dir = 0)
#' result$pval
circ_vtest <- function(alpha, dir, w = NULL, d = 0) {
  alpha <- as.numeric(alpha)

  if (is.null(w)) {
    w <- rep(1, length(alpha))
  } else {
    w <- as.numeric(w)
    if (length(w) != length(alpha)) {
      stop("Input dimensions do not match.")
    }
  }

  # Compute ingredients
  r <- circ_r(alpha, w, d)
  mu <- circ_mean(alpha, w)
  n <- sum(w)

  # Rayleigh's R (equ. 27.1)
  R <- n * r

  # V statistic (equ. 27.5)
  v <- R * cos(mu - dir)

  # u statistic (equ. 27.6)
  u <- v * sqrt(2 / n)

  # p-value from one-tailed normal approximation
  pval <- 1 - stats::pnorm(u)

  list(pval = pval, v = v)
}

#' Rayleigh Test for Non-uniformity
#'
#' Tests whether a sample of circular data is uniformly distributed.
#' This is a wrapper around the resultant length statistic.
#'
#' @param alpha numeric vector of angles in radians.
#' @param w optional weights (same length as \code{alpha}).
#' @return list with:
#'   \describe{
#'     \item{pval}{p-value of the Rayleigh test}
#'     \item{z}{Rayleigh's Z statistic}
#'   }
#' @references
#'   Zar JH (2010). Biostatistical Analysis. 5th ed. Prentice Hall.
#' @export
#' @examples
#' # Uniform data - expect high p-value
#' uniform_angles <- runif(100, 0, 2*pi)
#' circ_rayleigh(uniform_angles)
#'
#' # Concentrated data - expect low p-value
#' concentrated <- rnorm(100, mean = 0, sd = 0.3)
#' circ_rayleigh(concentrated)
circ_rayleigh <- function(alpha, w = NULL) {
  alpha <- as.numeric(alpha)
  n <- length(alpha)
  if (is.null(w)) {
    w <- rep(1, n)
  } else {
    w <- as.numeric(w)
    if (length(w) != n) stop("Input dimensions do not match.")
  }

  valid <- is.finite(alpha) & is.finite(w) & w > 0
  if (!any(valid)) return(list(pval = NA_real_, z = NA_real_))
  alpha <- alpha[valid]
  w <- w[valid]

  n_eff <- sum(w)
  if (!is.finite(n_eff) || n_eff <= 0) return(list(pval = NA_real_, z = NA_real_))

  r <- circ_r(alpha, w)
  R <- n_eff * r

  # Rayleigh's Z
  z <- R^2 / n_eff

  # p-value approximation (Zar, 2010)
  pval <- exp(-z) * (1 + (2 * z - z^2) / (4 * n_eff) -
                       (24 * z - 132 * z^2 + 76 * z^3 - 9 * z^4) / (288 * n_eff^2))
  pval <- max(min(pval, 1), 0)

  list(pval = pval, z = z)
}

#' FDR Correction (Benjamini-Hochberg)
#'
#' Wrapper around \code{stats::p.adjust} for FDR correction using the
#' Benjamini-Hochberg procedure. Provides interface similar to MATLAB \code{fdr_bh}.
#'
#' @param pvals numeric vector of p-values.
#' @param q desired false discovery rate (default 0.05).
#' @param method adjustment method passed to \code{p.adjust}. Default "BH".
#' @return list with:
#'   \describe{
#'     \item{h}{logical vector indicating significant tests after FDR correction}
#'     \item{crit_p}{critical p-value threshold (largest p-value deemed significant)}
#'     \item{adj_p}{adjusted p-values}
#'   }
#' @references
#'   Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate:
#'   a practical and powerful approach to multiple testing.
#'   Journal of the Royal Statistical Society B, 57, 289-300.
#' @export
#' @examples
#' pvals <- c(0.001, 0.01, 0.03, 0.04, 0.05, 0.1, 0.5)
#' result <- fdr_bh(pvals, q = 0.05)
#' result$h
fdr_bh <- function(pvals, q = 0.05, method = "BH") {
  adj_p <- stats::p.adjust(pvals, method = method)
  h <- adj_p <= q

  # Find critical p-value (largest original p deemed significant)
  if (any(h)) {
    crit_p <- max(pvals[h])
  } else {
    crit_p <- 0
  }

  list(h = h, crit_p = crit_p, adj_p = adj_p)
}
