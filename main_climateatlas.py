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
import os, psutil
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

files = fs.glob('resiclim/daily/total_precipitation*')

# define the output data path
outpath = '/projappl/project_2005030/climateatlas/out/'


## define the years
years = [2018, 2019, 2020]
years = np.arange(1991, 2021)


# select which indices to calculate
variables = {'growing_season_length':  True,
             'growing_degree_days':    True, 
             'freezing_degree_days':   True,
             'rain_on_snow':           False,
             'winter_warming':         True,
             'frost_growing_season':   True,
             'vapor_pressure_deficit': True,
             'heatwave_magnitude':     True, 
             'snow_season_length':     True}

# allocate dataarray dictionary
da_lists = dict((k, []) for k, v in variables.items() if v)

# allocate dataset dictionary
ds_out = dict.fromkeys(da_lists)

# variable short names
shortnames = {'growing_season_length':  'GSL',
              'growing_degree_days':    'GDD', 
              'freezing_degree_days':   'FDD',
              'rain_on_snow':           'ROS',
              'winter_warming':         'WW',
              'frost_growing_season':   'FGS',
              'vapor_pressure_deficit': 'VPD',
              'heatwave_magnitude':     'HWM', 
              'snow_season_length':     'SSL'}

# read heatwave threshold climatology
threshold_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/heatwave_threshold.nc')
T90p = threshold_ds.p90 - 273.15
    
# read 25th and 75th percentiles of annual maximum temperatures
percentiles_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/heatwave_25_75.nc')
p75max = percentiles_ds.p75_max - 273.15
p25max = percentiles_ds.p25_max - 273.15


for year in years:
    
    #### READ RAW INPUT ERA5-Land data-arrays #######
    
    da_t2mean_summer = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_temperature','DMEA')
    da_t2mean_winter = io_utils.read_daily_data_from_allas(fs, [year], 'winter', '2m_temperature','DMEA')
    da_t2max = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_temperature','DMAX')
    da_d2mean = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_dewpoint_temperature','DMEA')
    da_snowc = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'snow_cover','DMEA')
    # da_tp = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'total_precipitation','DSUM')
    da_sf = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'snowfall','DSUM')
    da_skt = io_utils.read_daily_data_from_allas(fs, [year], 'summer', 'skin_temperature','DMEA')


# da_t2mean = io_utils.read_test_data_from_lustre('2m_temperature', 'mean')
# da_t2max = io_utils.read_test_data_from_lustre('2m_temperature', 'max')
# da_d2mean = io_utils.read_test_data_from_lustre('2m_dewpoint_temperature', 'mean')
# da_tp = io_utils.read_test_data_from_lustre('total_precipitation',  'sum')
# da_sf = io_utils.read_test_data_from_lustre('snowfall',  'sum')
# da_snowc = io_utils.read_test_data_from_lustre('snow_cover',  'mean')
 

    # define the basevalue (threshold) in Celsius for growing season   
    basevalue = 5
    # define the rain threshold for ROS events in meters
    rain_threshold = 0.005
  
######### CALCULATE THE VARIOUS INDICES/VARIABLES

    # growing season length
    if variables["growing_season_length"]:
        gsl_tmp = indices.thermal_growing_season_length(da_t2mean_summer, basevalue)
        da_lists["growing_season_length"].append(gsl_tmp.compute())

    # growing degree days
    if variables["growing_degree_days"]:
        gdd_tmp = indices.thermal_growing_degree_days(da_t2mean_summer, basevalue)
        da_lists["growing_degree_days"].append(gdd_tmp.compute())
    
    # freezing degree days
    if variables["freezing_degree_days"]:
        fdd_tmp = indices.freezing_degree_days(da_t2mean_winter)
        da_lists["freezing_degree_days"].append(fdd_tmp.compute())

    # rain-on-snow events 
    if variables["rain_on_snow"]:
        ros_tmp = indices.rain_on_snow(da_tp, da_sf, da_snowc, rain_threshold)
        da_lists["rain_on_snow"].append(ros_tmp.compute())

    # # winter warming events 
    if variables["winter_warming"]:
        wwe_tmp = indices.winter_warming(da_t2mean_winter, da_snowc)
        da_lists["winter_warming"].append(wwe_tmp.compute())

    # # exposure to frost during growing season
    if variables["frost_growing_season"]:
        fgs_tmp = indices.frost_during_growing_season(da_t2mean_summer, da_skt, basevalue)
        da_lists["frost_growing_season"].append(fgs_tmp.compute())

    # # Vapour pressure deficit
    if variables["vapor_pressure_deficit"]:
        vpd_tmp = indices.vapour_pressure_deficit(da_t2mean_summer, da_d2mean)
        da_lists["vapor_pressure_deficit"].append(vpd_tmp.compute())

    # # Heatwave magnitude index
    if variables["heatwave_magnitude"]:
        hwi_tmp = indices.heatwave_magnitude_index(da_t2max, T90p, p75max, p25max)
        da_lists["heatwave_magnitude"].append(hwi_tmp.compute())
    
    # # Snow season length
    if variables["snow_season_length"]:
        ssl_tmp = indices.snow_season_length(da_snowc, )
        da_lists["snow_season_length"].append(ssl_tmp.compute())
    
    # print memory usage
    print('Memory usage:')
    print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 'MB')

print('All incides computed!')

# concatenate dataarrays
for var in da_lists:
    ds_out[var] = xr.concat(da_lists[var], dim='time').compute().to_dataset(name=shortnames[var])
    
    # add time attribute
    ds_out[var].time.attrs['long_name'] = "time"

    # add global attributes
    ds_out[var].attrs['Conventions'] = 'CF-1.7'
    ds_out[var].attrs['title'] = 'Bioclimatic indices'
    ds_out[var].attrs['Institution'] = 'Finnish Meteorological Institute'
    ds_out[var].attrs['source'] = 'ERA5-Land'
    ds_out[var].attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'
    
    # define outfile
    outfile = outpath  + 'resiclim_' + shortnames[var] + '.nc'
    # save the data as a netcdf file
    ds_out[var].to_netcdf(outfile, format='NETCDF4', encoding={'time': {'dtype': 'i4'}})
    
    

print('Done!')
