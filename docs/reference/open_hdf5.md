# Load HDF5 files (wrapper)

Thin wrapper around hdf5r::H5File to keep dependency optional.

## Usage

``` r
open_hdf5(path, mode = "r")
```

## Arguments

- path:

  file path to HDF5.

- mode:

  file mode (default "r").

## Value

H5File object.
