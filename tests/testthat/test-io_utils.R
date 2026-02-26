test_that("read_mat fails gracefully when dependency missing", {
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "R.matlab")) return(FALSE)
      base::requireNamespace(package, ...)
    },
    .package = "base"
  )
  expect_error(read_mat("dummy.mat"), "R.matlab required")
})

test_that("open_hdf5 fails gracefully when dependency missing", {
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "hdf5r")) return(FALSE)
      base::requireNamespace(package, ...)
    },
    .package = "base"
  )
  expect_error(open_hdf5("dummy.h5"), "hdf5r required")
})

# --- read_mat with package available ---

test_that("read_mat errors on non-existent file when R.matlab available", {
  skip_if_not_installed("R.matlab")
  suppressWarnings(expect_error(read_mat("/nonexistent/path/dummy.mat")))
})

# --- open_hdf5 with package available ---

test_that("open_hdf5 errors on non-existent file when hdf5r available", {
  skip_if_not_installed("hdf5r")
  expect_error(open_hdf5("/nonexistent/path/dummy.h5"))
})
