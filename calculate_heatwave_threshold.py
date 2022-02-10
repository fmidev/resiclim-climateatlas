#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 08:46:25 2022

This script calculates the 90th percentile threshold of climatological daily 
maximum temperatures, centred on a 31 d window, for the base period 1981–2010.
This threshold is used in the calculation of heatwave magnitude index.
The threshold is calculated for each grid point from the 30-year climatology
and 31-day moving window (31 x 30 values).

The heatwave magnitude index is based on this paper:
https://iopscience.iop.org/article/10.1088/1748-9326/ab6398/meta 

### NOTE!! ###
Currently the script calculates also 25th and 75th percentiles. This is probably
not necessary. 

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

def save_as_dataarray(array, doy_values, da_t2max):
    
    da = xr.DataArray(data=array, dims=["doy","latitude","longitude",],
                                           coords={'doy':doy_values,
                                                   'longitude': da_t2max.longitude, 
                                                   'latitude':da_t2max.latitude})
    return da


# define the output data path to which the resulting percentiles are saved
outpath = '/scratch/project_2005030/test/heatwave_25_75.nc'

## define the climatology years. 1981-2010 are used in the literature
years = np.arange(1991,2011)

# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 

fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})
files = fs.glob('resiclim/daily/2m_temperature_DMEA*')


#### READ RAW INPUT ERA5-Land #######
da_t2max = io_utils.read_daily_data_from_allas(fs, years, 'summer', '2m_temperature','DMAX')

# Width of the selection window (days)
struct = np.ones(31) 
    
# Arrays for saving the result
p90_array = np.empty((366, np.shape(da_t2max)[1], np.shape(da_t2max)[2]))
p75_array = np.empty((366, np.shape(da_t2max)[1], np.shape(da_t2max)[2]))
p25_array = np.empty((366, np.shape(da_t2max)[1], np.shape(da_t2max)[2]))

doy_values = np.unique(da_t2max['time.dayofyear'].values)
    
for day in doy_values:
    start2 = time.time()
    print('Calculating thresholds for day '+str(day))

    dayofyear = da_t2max['time.dayofyear'] == day
    
    # Pick 31 days around the target day from each year
    selection = binary_dilation(dayofyear, structure=struct)
    
    # ignore runtime warning
    warnings.simplefilter("ignore", category=RuntimeWarning)

    # Calculate the percentiles and save into arrays
    p90_array[int(day)-1,:, :] = np.nanpercentile(da_t2max.sel(time=selection), 90, axis=0)
    p75_array[int(day)-1,:, :] = np.nanpercentile(da_t2max.sel(time=selection), 75, axis=0)
    p25_array[int(day)-1,:, :] = np.nanpercentile(da_t2max.sel(time=selection), 25, axis=0)
    
    end2 = time.time()
    print('Day took '+str(np.round(end2 - start2, 2))+' seconds')


p90_da = save_as_dataarray(p90_array, doy_values, da_t2max)
p75_da = save_as_dataarray(p75_array, doy_values, da_t2max)
p25_da = save_as_dataarray(p25_array, doy_values, da_t2max)

### Create dataset for the output variables
ds_out = p90_da.to_dataset(name='p90')
ds_out['p75'] = p75_da
ds_out['p25'] = p25_da

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
ds_out.to_netcdf(outpath, format='NETCDF4')

print('Done!')

end = time.time()
print('Calculation took '+str(np.round(end - start, 2))+' seconds')