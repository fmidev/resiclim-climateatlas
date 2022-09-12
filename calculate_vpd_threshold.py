#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 08:46:25 2022

This script calculates the 90th percentile threshold of climatological daily 
VPD, centred on a 31 d window, for the base period 1981–2010.
This threshold is used in the calculation of VPD magnitude index.
The threshold is calculated for each grid point from the 30-year climatology
and 31-day moving window (31 x 30 values).

The calculation is done similarly as the heatwave threshold for HWMI.

@author: mprantan
"""

import os
import s3fs
import xarray as xr
import numpy as np
from datetime import datetime
import io_utils
from scipy.ndimage.morphology import binary_dilation
import warnings
import time


start = time.time()

def save_as_dataarray(array, doy_values, vpd):
    
    da = xr.DataArray(data=array, dims=["doy","latitude","longitude",],
                                           coords={'doy':doy_values,
                                                   'longitude': vpd.longitude, 
                                                   'latitude':vpd.latitude})
    return da


# define the output data path to which the resulting percentiles are saved
outpath = '/projappl/project_2005030/climateatlas/vpd_threshold.nc'

## define the climatology years. 1981-2010 are used in the literature
years = np.arange(1981,2011)

# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 

fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})
files = fs.glob('resiclim/daily/2m_temperature_DMEA*')


#### READ RAW INPUT ERA5-Land #######
da_t2mean = io_utils.read_daily_data_from_allas(fs, years, 'summer', '2m_temperature','DMEA')
da_d2mean = io_utils.read_daily_data_from_allas(fs, years, 'summer', '2m_dewpoint_temperature','DMEA')

# Select annual temperature
da_2t_annual = da_t2mean - 273.15
da_2d_annual = da_d2mean - 273.15

# Calculate Saturated Vapour Pressure in Pa using improved Magnus formula
VPsat = 610.94 * np.exp((17.625*da_2t_annual)/(da_2t_annual + 243.04))
    
# Calculate actual Vapour Pressure in Pa
VPair = 610.94 * np.exp((17.625*da_2d_annual)/(da_2d_annual + 243.04))
        
# Calculate the deficit
vpd = VPsat - VPair

# Width of the selection window (days)
struct = np.ones(31) 
    
# Array for saving the result
p90_array = np.empty((366, np.shape(vpd)[1], np.shape(vpd)[2]))


doy_values = np.unique(da_t2mean['time.dayofyear'].values)
    
for day in doy_values:
    start2 = time.time()
    print('Calculating thresholds for doy '+str(day))

    dayofyear = vpd['time.dayofyear'] == day
    
    # Pick 31 days around the target day from each year
    selection = binary_dilation(dayofyear, structure=struct)
    
    # ignore runtime warning
    warnings.simplefilter("ignore", category=RuntimeWarning)

    # Calculate the percentiles and save into arrays
    p90_array[int(day)-1,:, :] = np.nanpercentile(vpd.sel(time=selection), 90, axis=0)
    
    end2 = time.time()
    print('Calculation lasted '+str(np.round(end2 - start2, 2))+' seconds')


p90_da = save_as_dataarray(p90_array, doy_values, vpd)

### Create dataset for the output variables
ds_out = p90_da.to_dataset(name='p90')
ds_out = ds_out.compute()

# add time attribute
ds_out.doy.attrs['long_name'] = "Day of year"

# add global attributes
ds_out.attrs['Conventions'] = 'CF-1.7'
ds_out.attrs['title'] = 'Heatwave thresholds'
ds_out.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_out.attrs['source'] = 'ERA5-Land'
ds_out.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'


# save the data as a netcdf file
ds_out.to_netcdf(outpath, format='NETCDF4', encoding={'doy': {'dtype': 'i4'}})

print('Done!')

end = time.time()
print('Calculation took '+str(np.round(end - start, 2))+' seconds')