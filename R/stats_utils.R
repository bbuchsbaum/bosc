#' Pairwise Phase Consistency (1D)
#'
#' Bias-free PPC for angle data along a chosen dimension.
#'
#' @param data array of angles (radians).
#' @param dim integer; dimension along which to compute PPC.
#' @return array of PPC values with the chosen dimension collapsed.
#' @export
pairwise_phase_consistency <- function(data, dim = 1) {
  if (length(dim) != 1) stop("Only one dimension supported.")
  dims <- dim(data)
  if (length(dims) == 0) stop("Data must be at least 1D.")
  perm <- c(dim, setdiff(seq_along(dims), dim))
  dp <- aperm(data, perm)
  n <- dim(dp)[1]
  rest_dims <- dim(dp)
  rest_dims[1] <- 1
  N <- array(0, rest_dims)
  dumsum <- array(0, rest_dims)
  if (n < 2) return(array(NA_real_, rest_dims))
  for (p1 in 1:(n - 1)) {
    for (p2 in (p1 + 1):n) {
      dum <- cos(dp[p1, , drop = FALSE] - dp[p2, , drop = FALSE])
      idx <- !is.na(dum)
      N[idx] <- N[idx] + 1
      dum[is.na(dum)] <- 0
      dumsum <- dumsum + dum
    }
  }
  res <- dumsum / N
  aperm(res, order(perm))
}

#' Mann-Whitney U Z-score per element vs. reference repeats
#'
#' @param dat array of data.
#' @param ref array where first dimension is repetitions.
#' @param alpha significance threshold.
#' @return list with \code{Z} (same shape as \code{dat}) and \code{sgnf} (significance mask).
#' @export
u_score_matrix <- function(dat, ref, alpha = 0.05) {
  dshape <- dim(dat)
  if (is.null(dshape)) dshape <- length(dat)
  rshape <- dim(ref)
  if (is.null(rshape) || length(rshape) < 1) stop("ref must have first dim as repetitions.")
  n2 <- rshape[1]
  dat_vec <- as.vector(dat)
  ref_mat <- matrix(ref, nrow = n2)
  if (ncol(ref_mat) != length(dat_vec)) stop("ref dimensions must match dat beyond first dim.")
  # compare each column to corresponding dat element
  U <- colSums(ref_mat < matrix(dat_vec, nrow = n2, ncol = length(dat_vec), byrow = TRUE)) + 1
  mU <- n2 / 2
  sigU <- sqrt(n2 * (n2 + 2) / 12)
  Z <- (U - mU) / sigU
  Z_arr <- array(Z, dim = dshape)
  sgnf <- sign(Z_arr) * (abs(Z_arr) >= stats::qnorm(1 - alpha, 0, 1))
  list(Z = Z_arr, sgnf = sgnf)
}

#' Paired t-score along a dimension
#'
#' @param dat1 array.
#' @param dat2 optional array of same shape; if missing, one-sample t on \code{dat1}.
#' @param dim dimension to operate over.
#' @return list with \code{t} scores and \code{nu} degrees of freedom.
#' @export
paired_tscore <- function(dat1, dat2 = NULL, dim = 1) {
  if (!is.null(dat2) && !all(dim(dat1) == dim(dat2))) stop("dat1 and dat2 must match.")
  dat <- if (is.null(dat2)) dat1 else dat1 - dat2
  n <- dim(dat)[dim]
  tvals <- apply(dat, setdiff(seq_along(dim(dat)), dim), function(x) {
    xvec <- as.numeric(x)
    mn <- mean(xvec, na.rm = TRUE)
    sdv <- stats::sd(xvec, na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) return(NA_real_)
    mn / (sdv / sqrt(n))
  })
  # reshape back
  out_dim <- dim(dat)
  out_dim[dim] <- 1
  t_arr <- array(tvals, dim = out_dim)
  nu <- (n - 1) * array(1, dim = out_dim)
  list(t = t_arr, nu = nu)
}

#' Non-parametric p-values vs reference
#'
#' @param data numeric vector.
#' @param ref numeric vector of reference samples.
#' @param alpha significance level.
#' @return list with \code{h} significance (0/1) and \code{p} p-values.
#' @export
nonparam_pval <- function(data, ref, alpha = 0.05) {
  vals <- sort(ref)
  p <- numeric(length(data))
  h <- numeric(length(data))
  for (i in seq_along(data)) {
    mu <- stats::median(vals)
    if (data[i] < mu) {
      first_ge <- which(vals >= data[i])[1]
      dum <- if (is.na(first_ge)) length(vals) else max(0, first_ge - 1)
    } else {
      last_le <- max(which(vals <= data[i]))
      dum <- if (is.finite(last_le)) length(vals) - last_le else length(vals)
      dum <- max(0, dum)
    }
    p[i] <- dum / length(vals)
    h[i] <- as.numeric(p[i] <= alpha)
  }
  list(h = h, p = p)
}
