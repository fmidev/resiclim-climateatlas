#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 09:30:31 2022

This script calculates 25th and 75th percentiles of annual Tmax data from 1981-2010.
These percentiles are used in the calculation of heatwave magnitude index.
The percentiles are calculated for each grid point from the 30-year climatology (30 values).

The heatwave magnitude index is based on this paper:
https://iopscience.iop.org/article/10.1088/1748-9326/ab6398/meta 

@author: mprantan
"""

import os
import s3fs
import numpy as np
from datetime import datetime
import io_utils
import time



start = time.time()

# define the output path to which the resulting percentiles are saved
outpath = '/projappl/project_2005030/climateatlas/heatwave_25_75.nc'

## define the climatology years. 1981-2010 are used in the literature
years = np.arange(1981,2011)

# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 
fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})


#### READ RAW INPUT ERA5-Land #######
da_t2max = io_utils.read_daily_data_from_allas(fs, years, 'summer', '2m_temperature','DMAX')

# calculate annual maxima
print('Calculating annual maximum values...')
t2max_annual = da_t2max.groupby(da_t2max.time.dt.year).max().compute()

# calculate the percentiles 
print('Calculating 25th and 75th percentiles...')
p75_max = t2max_annual.quantile(0.75, dim='year')
p25_max = t2max_annual.quantile(0.25, dim='year')


### Create dataset for the output variables
ds_out = p75_max.to_dataset(name='p75_max')
ds_out['p25_max'] = p25_max

# compute the variables
ds_out = ds_out.compute()


# add global attributes
ds_out.attrs['Conventions'] = 'CF-1.7'
ds_out.attrs['title'] = '25th and 75th percentiles of annual maximum temperatures'
ds_out.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_out.attrs['source'] = 'ERA5-Land'
ds_out.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'


# save the data as a netcdf file
ds_out.to_netcdf(outpath, format='NETCDF4')

print('Done!')

end = time.time()
print('Calculation took '+str(np.round(end - start, 2))+' seconds')