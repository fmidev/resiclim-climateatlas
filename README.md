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
python=3.7.6
numpy=1.18.1
scipy=1.4.1
pandas=0.25.3
geopandas=0.9.0
fiona=1.8.13
gdal=3.0.4
matplotlib=3.1.2
cartopy=0.17.0
seaborn=0.11.1
scikit-learn=0.21.3
xgboost=1.4.2
xarray=0.13.0
bottleneck=1.3.0
dask=2.8.0
s3fs=0.4.2
boto3=1.19.8
netCDF4=1.5.3
h5netcdf=0.11.0
ecmwf-api-client=1.6.1
cdsapi=0.5.1
spyder=5.2.1
