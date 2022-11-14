# ARCLIM climate atlas
Climate and event-type indices derived and published in the RESICLIM project

## Python files

### Main program
The ARCLIM variables are calculated in `main_climateatlas.py`.
The code reads daily ERA5-Land data from CSC Allas object storage system and calculates
the ARCLIM indices. The ARCLIM indices are saved in both NETCDF4 and GeoTIFF format.

### Calculation of the indices
The module to calculate the various indices can be found from `indices.py`.

### Input/output
IO-utils are written in `io_utils.py`.

### Downloading ERA5-Land to Allas
These codes are in `allas_utils.py`.

## Requirements for python packages
The yml file consisting of the python packages is provided in the repository (resiclim.yml)
