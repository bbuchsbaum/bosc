test_that("read_mat fails gracefully when dependency missing", {
  if (requireNamespace("R.matlab", quietly = TRUE)) skip("R.matlab installed")
  expect_error(read_mat("dummy.mat"), "R.matlab required")
})

test_that("open_hdf5 fails gracefully when dependency missing", {
  if (requireNamespace("hdf5r", quietly = TRUE)) skip("hdf5r installed")
  expect_error(open_hdf5("dummy.h5"), "hdf5r required")
})
