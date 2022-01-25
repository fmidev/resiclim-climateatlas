#!/bin/sh -l
# 
#
#
# Commands to setup and istall the Miniconda environment and to 
# initialize the S3 storage system in CSC Allas for the RESICLIM project
#
# Some of these commands work only in the CSC Puhti supercomputer
#
# This file must be run only once in the beginning of the project
#
#
#  S3 Endpoint: a3s.fi
#  DNS-style bucket+hostname:port template for accessing a bucket: %(bucket)s.a3s.fi
#  website_endpoint = http://%(bucket)s.s3-website-%(location)s.amazonaws.com/
#
#  The public address of data (if exists) will follow this format:
#  https://resiclim.a3s.fi/some_folder_path/data.nc



cd /projappl/project_2005404

# Fetch the installation package
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh 

# When executing the installer, change the default installation path to /projappl/project_2005404/miniconda
bash Miniconda3-py37_4.8.2-Linux-x86_64.sh 


export PATH="/projappl/project_2005404/miniconda/bin:$PATH"
export LD_LIBRARY_PATH="/projappl/project_2005404/miniconda/lib"


# Install additional packages
conda install -c conda-forge \
        python=3.7.6 \
        numpy=1.18.1 \
        scipy=1.4.1 \
        pandas=0.25.3 \
        geopandas=0.9.0 \
        fiona=1.8.13 \
        gdal=3.0.4 \
        matplotlib=3.1.2 \
        cartopy=0.17.0 \
        seaborn=0.11.1 \
        scikit-learn=0.21.3 \
        xgboost=1.4.2 \
        xarray=0.13.0 \
        zarr=2.10.3 \
        bottleneck=1.3.0 \
        dask=2.8.0 \
        s3fs=0.4.2 \
        boto3=1.19.8 \
        netCDF4=1.5.3 \
        h5netcdf=0.11.0 \
        ecmwf-api-client=1.6.1 \
        cdsapi=0.5.1 \
        spyder=5.2.1






mkdir -p /users/kamarain/resiclim
cd /users/kamarain/resiclim

# Allas environment
module load allas

# Create the configuration file for the s3cmd tool and rename it
allas-conf project_2005404 --mode s3cmd
mv $HOME/.s3cfg $HOME/.s3cfg_resiclim


# Use RESICLIM credentials 
alias s3cmd_resic='s3cmd --config=/users/kamarain/.s3cfg_resiclim'





# Create a bucket for storing the ERA5 files 
s3cmd_resic mb --acl-public   s3://resiclim


# s3cmd_resic put 2m_temperature_era5Land_202001.nc s3://resiclim-era5/2m_temperature/ -P
# https://resiclim.a3s.fi/2m_temperature/2m_temperature_era5Land_202001.nc



# Delete everything
#s3cmd_resic rm  s3://resiclim/2m_temperature/*
#s3cmd_resic rm  s3://resiclim/total_precipitation/*
#s3cmd_resic rb  s3://resiclim


