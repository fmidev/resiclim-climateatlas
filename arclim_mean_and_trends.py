#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 14:02:29 2022

This script calculates mean conditions 1991-2020 and trends for 1951-2021
from ARCLIM annual values. The trends are calculated onwards from 1951 
because the year 1951 is the first complete year in winter variables.

The trends are calculated using Theil slope method.


@author: mprantan
"""

import os, glob
import xarray as xr
import numpy as np
from datetime import datetime
import pymannkendall as mk
import io_utils
import rioxarray

# define the input path
path_for_annual_data = '/projappl/project_2005030/climateatlas/out/annual/nc/'

# define the output data path
outpath = '/projappl/project_2005030/climateatlas/out/'

# define files
variables_in_order = ['GSL','GDD','FGS','FDD','ROS','WWE','WWI','HWMI','VPDI',
                      'SWI','SSL','SSO','SSE','HWE','TAVG','PRA','SFA','WSA']

# append the path to the variable names
files = [path_for_annual_data + "arclim_"+ s +".nc" for s in variables_in_order]


# Allocate lists for the variables
list_of_trends = []
list_of_means = []
list_of_pvals = []

# Print every 10000 lines.
LOG_EVERY_N = 10000

# loop over all files
for i, f in enumerate(files):
    
    # open dataarray and select 1951-2021 period
    da = xr.open_dataarray(f).sel(time=slice('1951','2021'))
    
     # 1/np.nan field (land sea mask)
    ls_mask = da.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
    
    # define the variable name
    var = da.name
    
    # years
    years = da.time.values
    year1= str(years[0])
    year2= str(years[-1])
    
    # Reshape to an array with as many rows as years and as many columns as there are pixels
    vals = da.values.reshape(len(years), -1)
    slopes = np.zeros(( np.shape(vals)[1]))
    slopes[:] = np.nan
    p_vals = np.zeros(( np.shape(vals)[1]))
    p_vals[:] = np.nan

    size_of_array = np.shape(vals)[1]
    
    #loop over all grid cells and calculate slope if the cell is land
    for p in np.arange(0, np.shape(vals)[1]):
        if (p % LOG_EVERY_N) == 0:
            print(i, var, np.round(p/size_of_array * 100,1), '%')
        y = vals[:, p]
        # check if the grid cell values are not nans
        if all(np.isfinite(y)):
            # Perform Mann-Kendall test
            res = mk.original_test(y)
            # Slope of the trend
            slopes[p] = res.slope
            # p-value of the trend
            p_vals[p] = res.p
    
    # reshape back to 2d arrays
    trends = slopes.reshape(da.values.shape[1], da.values.shape[2])
    p_values = p_vals.reshape(da.values.shape[1], da.values.shape[2])
     
    # make labelled xarray and multiply by the land-sea mask       
    trends_da = da.sel(time=2021).copy(data=trends) * ls_mask
    p_values_da = da.sel(time=2021).copy(data=p_values) * ls_mask
   
    # add attributes
    trends_da.attrs['long_name'] = da.long_name + ' '+year1+'-'+year2+' trend'
    trends_da.attrs['short_name'] = da.short_name
    p_values_da.attrs['long_name'] = da.long_name + ' '+year1+'-'+year2+' p-value'
    p_values_da.attrs['short_name'] = da.short_name
    
    # save trend array
    list_of_trends.append(trends_da.rename(var))
    list_of_pvals.append(p_values_da.rename(var))
    
    # calculate mean conditions over 1991-2020
    mean_da = da.sel(time=slice('1991','2020')).mean(dim='time')
    mean_da.attrs['long_name'] = da.long_name + ' 1991-2020'
    mean_da.attrs['short_name'] = da.short_name
    
    # save mean array
    list_of_means.append(mean_da.rename(var))
    
  
# make datasets
ds_trend = list_of_trends[0].to_dataset()
ds_pvals = list_of_pvals[0].to_dataset()
ds_mean = list_of_means[0].to_dataset()



for d in list_of_trends[1:]:
    ds_trend[d.name] = d

    
# add global attributes
ds_trend.attrs['title'] = 'ARCLIM Bioclimatic indices'
ds_trend.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_trend.attrs['source'] = 'ERA5-Land'
ds_trend.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'
    
# define outfile
outfile = outpath  + 'arclim_trends.nc'
outfile_tiff = outpath  + 'arclim_trends.tif'
# save the data as a netcdf file
ds_trend.to_netcdf(outfile, format='NETCDF4')

# convert nc to geotiff and save
geotiff = io_utils.da_to_geotiff(ds_trend.to_array(dim='variable'))
geotiff.to_raster(outfile_tiff)




for d in list_of_pvals[1:]:
    ds_pvals[d.name] = d

    
# add global attributes
ds_pvals.attrs['title'] = 'ARCLIM Bioclimatic indices'
ds_pvals.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_pvals.attrs['source'] = 'ERA5-Land'
ds_pvals.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'
    
# define outfile
outfile = outpath  + 'arclim_pvalues.nc'
outfile_tiff = outpath  + 'arclim_pvalues.tif'

# save the data as a netcdf file
ds_pvals.to_netcdf(outfile, format='NETCDF4')

# convert nc to geotiff and save
geotiff = io_utils.da_to_geotiff(ds_trend.to_array(dim='variable'))
geotiff.to_raster(outfile_tiff)



for d in list_of_means[1:]:
    ds_mean[d.name] = d
    
# add global attributes
ds_mean.attrs['title'] = 'ARCLIM Bioclimatic indices'
ds_mean.attrs['Institution'] = 'Finnish Meteorological Institute'
ds_mean.attrs['source'] = 'ERA5-Land'
ds_mean.attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'
    
# define outfile
outfile = outpath  + 'arclim_means.nc'
outfile_tiff = outpath  + 'arclim_means.tif'
# save the data as a netcdf file
ds_mean.to_netcdf(outfile, format='NETCDF4')

# convert nc to geotiff and save
geotiff = io_utils.da_to_geotiff(ds_trend.to_array(dim='variable'))
geotiff.to_raster(outfile_tiff)
