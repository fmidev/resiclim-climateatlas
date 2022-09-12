#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 09:30:31 2022

This script calculates 25th and 75th percentiles of annual VPD data from 1981-2010.
These percentiles are used in the calculation of VPD magnitude index.
The percentiles are calculated for each grid point from the 30-year climatology (30 values).

VPD index is similar to heatwave magnitude index, which is in turn based on this paper:
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
outpath = '/projappl/project_2005030/climateatlas/vpd_25_75.nc'

## define the climatology years. 1981-2010 are used in the literature
years = np.arange(1981,2011)

# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 
fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})


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

# calculate annual maxima
print('Calculating annual maximum values...')
vpd_annual = vpd.groupby(vpd.time.dt.year).max().compute()

# calculate the percentiles 
print('Calculating 25th and 75th percentiles...')
p75_max = vpd_annual.quantile(0.75, dim='year').drop_vars('quantile')
p25_max = vpd_annual.quantile(0.25, dim='year').drop_vars('quantile')


### Create dataset for the output variables
ds_out = p75_max.to_dataset(name='p75_max')
ds_out['p25_max'] = p25_max

# compute the variables
ds_out = ds_out.compute()


# add global attributes
ds_out.attrs['Conventions'] = 'CF-1.7'
ds_out.attrs['title'] = '25th and 75th percentiles of annual VPD maxima'
ds_out.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_out.attrs['source'] = 'ERA5-Land'
ds_out.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'


# save the data as a netcdf file
ds_out.to_netcdf(outpath, format='NETCDF4')

print('Done!')

end = time.time()
print('Calculation took '+str(np.round(end - start, 2))+' seconds')