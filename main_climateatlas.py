#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 14:02:29 2022

This is the main script to calculate various bioclimatic indicators/variables
in RESICLIM project. 

The raw input ERA5-Land data is read from CSC Allas file system.

The indicators are calculated in a module called "indices".

Currently the output data are saved to one netcdf file

@author: mprantan
"""
import os
import s3fs
import xarray as xr
import numpy as np
from datetime import datetime
import indices
import io_utils


# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 

fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})

files = fs.glob('resiclim/daily/2m_temperature_DMEA*')

# define the output data path
outpath = '/scratch/project_2005030/test/test_indices.nc'


## define the years
years = [2018, 2019, 2020]
years = np.arange(2011, 2021)


gsl = []
gdd = []
ros = []
wwe = []
fgs = []
vpd = []
hwi = []

# read heatwave threshold climatology
threshold_ds = xr.open_dataset('/scratch/project_2005030/test/heatwave_threshold.nc')
T90p = threshold_ds.p90 - 273.15
    
# read 25th and 75th percentiles of annual maximum temperatures
percentiles_ds = xr.open_dataset('/scratch/project_2005030/test/heatwave_25_75.nc')
p75max = percentiles_ds.p75_max - 273.15
p25max = percentiles_ds.p25_max - 273.15


for year in years:
    
    #### READ RAW INPUT ERA5-Land data-arrays #######
    
    da_t2mean_summer = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_temperature','DMEA')
    da_t2mean_winter = io_utils.read_daily_data_from_allas(fs, [year], 'winter', '2m_temperature','DMEA')

    da_t2max = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_temperature','DMAX')
    da_d2mean = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_dewpoint_temperature','DMEA')
    da_snowc = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'snow_cover','DMEA')
    da_tp = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'total_precipitation','DMAX')
    da_sf = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'snowfall','DMAX')
    da_skt = io_utils.read_daily_data_from_allas(fs, [year], 'summer', 'skin_temperature','DMEA')


# da_t2mean = io_utils.read_test_data_from_lustre('2m_temperature', 'mean')
# da_t2max = io_utils.read_test_data_from_lustre('2m_temperature', 'max')
# da_d2mean = io_utils.read_test_data_from_lustre('2m_dewpoint_temperature', 'mean')
# da_tp = io_utils.read_test_data_from_lustre('total_precipitation',  'sum')
# da_sf = io_utils.read_test_data_from_lustre('snowfall',  'sum')
# da_snowc = io_utils.read_test_data_from_lustre('snow_cover',  'mean')
 

    # define the basevalue (threshold) in Celsius for growing season   
    basevalue = 3
    # define the rain threshold for ROS events in meters
    rain_threshold = 0.01
  
######### CALCULATE THE VARIOUS INDICES/VARIABLES

    # growing season length
    gsl_tmp = indices.thermal_growing_season_length(da_t2mean_summer, basevalue)
    gsl.append(gsl_tmp)

    # growing degree days 
    gdd_tmp = indices.thermal_growing_degree_days(da_t2mean_summer, basevalue)
    gdd.append(gdd_tmp)

    # rain-on-snowevents 
    ros_tmp = indices.rain_on_snow(da_tp, da_sf, da_snowc, rain_threshold)
    ros.append(ros_tmp)

    # winter warming events 
    wwe_tmp = indices.winter_warming(da_t2mean_winter, da_snowc)
    wwe.append(wwe_tmp)

    # exposure to frost during growing season
    fgs_tmp = indices.frost_during_growing_season(da_t2mean_summer, basevalue)
    fgs.append(fgs_tmp)

    # Vapour pressure deficit 
    vpd_tmp = indices.vapour_pressure_deficit(da_t2mean_summer, da_d2mean)
    vpd.append(vpd_tmp)

    # Heatwave magnitude index
    hwi_tmp = indices.heatwave_magnitude_index(da_t2max, T90p, p75max, p25max)
    hwi.append(hwi_tmp)

print('All incides calculated lazily')

# concatenate datasets
gsl = xr.concat(gsl, dim='time')
gdd = xr.concat(gdd, dim='time')
ros = xr.concat(ros, dim='time')
wwe = xr.concat(wwe, dim='time')
fgs = xr.concat(fgs, dim='time')
vpd = xr.concat(vpd, dim='time')
hwi = xr.concat(hwi, dim='time')


### Create dataset for the output variables
ds_out = gsl.to_dataset(name='gsl')
ds_out['gdd'] = gdd
ds_out['ros'] = ros
ds_out['wwe'] = wwe
ds_out['fgs'] = fgs
ds_out['vpd'] = vpd
ds_out['hwi'] = hwi


print('Computing incides...')
### Compute indices ###
ds_out = ds_out.compute()

print('All incides computed')


# add time attribute
ds_out.time.attrs['long_name'] = "time"

# add global attributes
ds_out.attrs['Conventions'] = 'CF-1.7'
ds_out.attrs['title'] = 'Bioclimatic indices'
ds_out.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_out.attrs['source'] = 'ERA5-Land'
ds_out.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'


# save the data as a netcdf file
ds_out.to_netcdf(outpath, format='NETCDF4', encoding={'time': {'dtype': 'i4'}})

print('Done!')
