#' Connected component labeling for thresholded data
#'
#' Thresholds 2D data and labels connected components (8-connectivity).
#'
#' @param data matrix of values.
#' @param thres numeric threshold; if NULL, assumes binary input.
#' @return integer matrix of labels (0 for below threshold).
#' @export
extract_clusters <- function(data, thres = NULL) {
  dat <- data
  if (!is.null(thres)) {
    dat <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    dat[data >= thres] <- 1
  }
  # pad boundaries with zeros
  dat[1, ] <- 0
  dat[nrow(dat), ] <- 0
  dat[, 1] <- 0
  dat[, ncol(dat)] <- 0

  label <- 2L
  labels <- matrix(0L, nrow = nrow(dat), ncol = ncol(dat))
  parent <- list()

  for (r in 2:(nrow(dat) - 1)) {
    for (c in 2:(ncol(dat) - 1)) {
      if (dat[r, c] != 0) {
        neighbors <- c(labels[r - 1, (c - 1):(c + 1)], labels[r, c - 1])
        neighbors <- neighbors[neighbors != 0]
        if (length(neighbors) == 0) {
          labels[r, c] <- label
          parent[[label]] <- label
          label <- label + 1L
        } else {
          labels[r, c] <- min(neighbors)
          ul <- unique(neighbors)
          for (u in ul) parent[[u]] <- min(c(parent[[u]], ul))
        }
      }
    }
  }
  # second pass: flatten parents
  if (length(parent) > 0) {
    for (r in 2:(nrow(dat) - 1)) {
      for (c in 2:(ncol(dat) - 1)) {
        if (labels[r, c] != 0) labels[r, c] <- parent[[labels[r, c]]]
      }
    }
  }
  labels
}

#' Detect clusters with optional splitting
#'
#' @param data matrix of primary values.
#' @param threshold numeric; binarization threshold.
#' @param method "extract" or "watershed".
#' @param smooth optional numeric or length-2 vector for Gaussian smoothing (pixels).
#' @param split_clusters logical; attempt to split multi-peak clusters.
#' @param peak_height numeric vector of relative peak heights for splitting.
#' @param min_new_cluster_size minimum fraction of original cluster to keep split.
#' @param freqweight,timeweight weights for distance reassignment.
#' @param freq_axis,time_axis vectors for frequencies/times corresponding to rows/cols.
#' @param data2 optional secondary matrix for peak value extraction.
#' @param dataRef optional reference matrix for computing reference statistics per cluster.
#' @param qval optional FDR q for p-value adjustments.
#' @return list of cluster structs with coordinates and stats.
#' @export
detect_clusters <- function(data,
                            threshold,
                            method = c("watershed", "extract"),
                            smooth = NULL,
                            split_clusters = FALSE,
                            peak_height = seq(0.05, 0.5, by = 0.05),
                            min_new_cluster_size = 0.1,
                            freqweight = 1,
                            timeweight = 1,
                            freq_axis = NULL,
                            time_axis = NULL,
                            data2 = NULL,
                            dataRef = NULL,
                            qval = NULL) {
  method <- match.arg(method)
  if (is.null(freq_axis)) freq_axis <- seq_len(nrow(data))
  if (is.null(time_axis)) time_axis <- seq_len(ncol(data))

  if (!is.null(smooth)) {
    if (!requireNamespace("imager", quietly = TRUE)) {
      stop("The 'imager' package is required for Gaussian smoothing. ",
           "Install it with install.packages('imager').", call. = FALSE)
    }
    sigma <- if (length(smooth) == 1) c(smooth, smooth) else smooth
    img <- imager::as.cimg(t(data))
    img <- imager::isoblur(img, sigma = sigma[1])
    sdata <- t(as.matrix(img))
  } else {
    sdata <- data
  }

  if (method == "extract") {
    clusDatPlus <- extract_clusters(sdata, threshold)
    clusDatMin <- extract_clusters(-sdata, threshold)
    clusDat <- -clusDatMin + clusDatPlus
  } else {
    # watershed via imager
    if (!requireNamespace("imager", quietly = TRUE)) {
      stop("The 'imager' package is required for the 'watershed' method. ",
           "Install it with install.packages('imager').", call. = FALSE)
    }
    posmask <- sdata
    posmask[sdata < threshold] <- 0
    posimg <- imager::as.cimg(t(posmask))
    clusPos <- t(as.matrix(imager::watershed(posimg)))

    negmask <- sdata
    negmask[-sdata > threshold] <- 0
    negimg <- imager::as.cimg(t(negmask))
    clusNeg <- -t(as.matrix(imager::watershed(negimg)))
    clusDat <- clusPos + clusNeg
    clusDat[!is.finite(clusDat)] <- 0
  }

  # optional splitting (simplified heuristic)
  if (split_clusters) {
    labels <- setdiff(unique(clusDat), 0)
    for (lab in labels) {
      mask <- clusDat == lab
      mm <- max(abs(data[mask]))
      # Use sign factor to preserve cluster polarity (matches MATLAB fact = 1 or -1)
      fact <- if (lab > 0) 1 else -1
      for (ph in peak_height) {
        dataRed <- fact * sdata - mm * ph
        dataRed[!mask] <- 0
        subClus <- extract_clusters(dataRed, thres = 0)
        newLabs <- setdiff(unique(subClus), 0)
        sizes <- vapply(newLabs, function(l) sum(subClus == l), numeric(1))
        if (sum(sizes / sum(mask) > min_new_cluster_size) >= 2) {
          keepLabs <- newLabs[sizes / sum(mask) > min_new_cluster_size]
          subMask <- subClus
          subMask[!subMask %in% keepLabs] <- 0
          # reassign remainder
          missing <- mask & subMask == 0
          if (any(missing)) {
            coords <- which(missing, arr.ind = TRUE)
            distmat <- matrix(Inf, nrow = nrow(coords), ncol = length(keepLabs))
            for (j in seq_along(keepLabs)) {
              pts <- which(subMask == keepLabs[j], arr.ind = TRUE)
              for (k in seq_len(nrow(coords))) {
                df <- (freq_axis[coords[k, 1]] - freq_axis[pts[, 1]]) / (max(freq_axis) - min(freq_axis))
                dt <- (time_axis[coords[k, 2]] - time_axis[pts[, 2]]) / (max(time_axis) - min(time_axis))
                distmat[k, j] <- min(sqrt(freqweight * df ^ 2 + timeweight * dt ^ 2))
              }
            }
            assign <- apply(distmat, 1, which.min)
            for (k in seq_len(nrow(coords))) {
              subMask[coords[k, 1], coords[k, 2]] <- keepLabs[assign[k]]
            }
          }
          offset <- if (lab > 0) max(clusDat) + 1 else min(clusDat) - 1
          clusDat[mask] <- offset + subMask[mask]
          break
        }
      }
    }
  }

  clusLabels <- setdiff(unique(clusDat), 0)
  if (length(clusLabels) == 0) return(list())

  out <- vector("list", length(clusLabels))
  for (i in seq_along(clusLabels)) {
    clusid <- clusLabels[i]
    locs <- which(clusDat == clusid, arr.ind = TRUE)
    peak_src <- if (!is.null(data2)) data2 else data
    peak_val <- peak_src[clusDat == clusid]
    peak_idx <- which.max(abs(peak_val))
    sumZ <- sum(data[clusDat == clusid])
    pAdj <- NA_real_
    if (!is.null(qval)) {
      p <- stats::pnorm(-abs(data))
      pAdj <- matrix(stats::p.adjust(p, method = "BH"), nrow = nrow(data))
    }
    out[[i]] <- list(
      blobID = clusid,
      times = time_axis[locs[, 2]],
      freqs = freq_axis[locs[, 1]],
      Zscores = data[clusDat == clusid],
      CoMtime = time_axis[round(mean(locs[, 2]))],
      CoMfreq = freq_axis[round(mean(locs[, 1]))],
      pAdj = if (is.matrix(pAdj)) pAdj[clusDat == clusid] else NA_real_,
      peakZscore = peak_val[peak_idx],
      sumZscore = sumZ,
      peakpAdj = if (is.matrix(pAdj)) pAdj[which(clusDat == clusid)[peak_idx]] else NA_real_
    )
    if (!is.null(dataRef)) {
      refvals <- dataRef[clusDat == clusid]
      out[[i]]$peakZscoreRef <- refvals[peak_idx]
      out[[i]]$maxZscoreRef <- max(refvals, na.rm = TRUE)
      out[[i]]$avgZscoreRef <- mean(refvals, na.rm = TRUE)
    }
  }
  out
}
