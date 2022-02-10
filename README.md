# resiclim-climateatlas
Climate indices derived and published in the RESICLIM project

## Python files

### Main program
The main analysis can be done with `main_climateatlas.py`.
The code is tailored to use CSC Allas file system.

### Calculation of the indices
The various bioclimatic indices are calculated with functions which can be found from `indices.py`.

### Input/output
IO-utils are written in `io_utils.py`.

### Downloading ERA5-Land to Allas
These codes are in `allas_utils.py`.

## Requirements
Python 3+, Xarray, Pandas, Scipy, ...
