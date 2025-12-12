#' Load MATLAB .mat files (wrapper)
#'
#' Thin wrapper around R.matlab::readMat to keep dependency optional.
#'
#' @param path file path to .mat file.
#' @return list from readMat.
#' @export
read_mat <- function(path) {
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("R.matlab required to read .mat files.")
  }
  R.matlab::readMat(path)
}

#' Load HDF5 files (wrapper)
#'
#' Thin wrapper around hdf5r::H5File to keep dependency optional.
#'
#' @param path file path to HDF5.
#' @param mode file mode (default "r").
#' @return H5File object.
#' @export
open_hdf5 <- function(path, mode = "r") {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r required to read HDF5 files.")
  }
  hdf5r::H5File$new(path, mode = mode)
}
