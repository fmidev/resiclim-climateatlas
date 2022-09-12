#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 14:02:29 2022

This is the main script to calculate various bioclimatic indicators/variables
for ARCLIM dataset in RESICLIM project. 

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
import rioxarray


# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 

fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})

files = fs.glob('resiclim/daily/*10m_v_component*DMEA*2010*')

# define the output data path
outpath = '/projappl/project_2005030/climateatlas/out/annual/'


# define the years
# years = [2015]
years = np.arange(1950, 2022)

# Select if you want to calculate all variables. If False, select the variables
# from below
calculate_all = True

# select which indices to calculate
variables = {'growing_season_length':    False,
             'growing_degree_days':      False, 
             'frost_growing_season':     False,
             'freezing_degree_days':     False,
             'rain_on_snow_freezing':    False,
             'rain_on_snow_thawing':     False,
             'winter_warming_events':    False,
             'winter_warming_intensity': False,
             'heatwave_magnitude':       True,
             'vpd_magnitude_index':      True,
             'summer_warmth_index':      True,
             'snow_season_length':       False,
             'high_wind_events':         False, 
             'annual_mean_temperature':  False,
             'annual_precipitation':     False,
             'annual_snowfall':          False,
             'average_wind_speed':       False,
             }

if calculate_all:
    variables = dict((k, True) for k, v in variables.items())

# allocate dataarray dictionary
da_lists = dict((k, []) for k, v in variables.items() if v)

# allocate dataset dictionary
ds_out = dict.fromkeys(da_lists)

# variable short names
shortnames = {'growing_season_length':    'GSL',
              'growing_degree_days':      'GDD', 
              'frost_growing_season':     'FGS',
              'freezing_degree_days':     'FDD',
              'rain_on_snow_freezing':    'ROSF',
              'rain_on_snow_thawing':     'ROST',
              'winter_warming_events':    'WWE',
              'winter_warming_intensity': 'WWI',
              'heatwave_magnitude':       'HWMI', 
              'vpd_magnitude_index':      'VPDI',
              'summer_warmth_index':      'SWI',
              'snow_season_length':       'SSL',
              'high_wind_events':         'HWE',
              'annual_mean_temperature':  'TAVG',
              'annual_precipitation':     'PRA',
              'annual_snowfall':          'SFA',
              'average_wind_speed':       'WSA',}

# read heatwave threshold climatology
if variables["heatwave_magnitude"]:
    threshold_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/heatwave_threshold.nc')
    T90p = threshold_ds.p90 - 273.15
    
    # read 25th and 75th percentiles of annual maximum temperatures
    percentiles_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/heatwave_25_75.nc')
    p75max = percentiles_ds.p75_max - 273.15
    p25max = percentiles_ds.p25_max - 273.15
    
 # read vpd threshold climatology
if variables["vpd_magnitude_index"]:
    threshold_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/vpd_threshold.nc')
    T90p_vpd = threshold_ds.p90
    
    # read 25th and 75th percentiles of annual maximum temperatures
    percentiles_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/vpd_25_75.nc')
    p75vpd = percentiles_ds.p75_max
    p25vpd = percentiles_ds.p25_max

# define the basevalue (threshold) in Celsius for growing season   
basevalue = 5
# define the rain threshold for ROS events in meters
rain_threshold = 0.005

# loop over the years
for year in years:
    
    #### READ PRE-PROCESSED INPUT ERA5-Land data-arrays #######
    
    da_t2mean_summer = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_temperature','DMEA')
    da_t2mean_winter = io_utils.read_daily_data_from_allas(fs, [year], 'winter', '2m_temperature','DMEA')
    da_t2max = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_temperature','DMAX')
    da_d2mean = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '2m_dewpoint_temperature','DMEA')
    da_snowc = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'snow_cover','DMEA')
    da_tp_summer = io_utils.read_daily_data_from_allas(fs, [year], 'summer', 'total_precipitation','DSUM')
    da_tp_winter = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'total_precipitation','DSUM')
    da_sf = io_utils.read_daily_data_from_allas(fs, [year], 'winter', 'snowfall','DSUM')
    da_skt = io_utils.read_daily_data_from_allas(fs, [year], 'summer', 'skin_temperature','DMIN')
    da_u10mean = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '10m_u_component_of_wind','DMEA')
    da_v10mean = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '10m_v_component_of_wind','DMEA')
    da_u10max = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '10m_u_component_of_wind','DMAX')
    da_v10max = io_utils.read_daily_data_from_allas(fs, [year], 'summer', '10m_v_component_of_wind','DMAX')
 

######### CALCULATE THE VARIOUS INDICES/VARIABLES

    # growing season length
    if variables["growing_season_length"]:
        gsl_tmp = indices.thermal_growing_season_length(da_t2mean_summer, basevalue)
        da_lists["growing_season_length"].append(gsl_tmp.compute())

    # growing degree days
    if variables["growing_degree_days"]:
        gdd_tmp = indices.thermal_growing_degree_days(da_t2mean_summer, basevalue)
        da_lists["growing_degree_days"].append(gdd_tmp.compute())
        
    # # exposure to frost during growing season
    if variables["frost_growing_season"]:
        fgs_tmp = indices.frost_during_growing_season(da_t2mean_summer, da_skt, basevalue)
        da_lists["frost_growing_season"].append(fgs_tmp.compute())
    
    # freezing degree days
    if variables["freezing_degree_days"]:
        fdd_tmp = indices.freezing_degree_days(da_t2mean_winter)
        da_lists["freezing_degree_days"].append(fdd_tmp.compute())
        
    # # Snow season length
    if variables["snow_season_length"]:
        ssl_tmp = indices.snow_season_length(da_snowc, )
        da_lists["snow_season_length"].append(ssl_tmp.compute())

    # rain-on-snow events 
    if variables["rain_on_snow_freezing"]:
        rosf_tmp, rost_tmp = indices.rain_on_snow(da_tp_winter, da_sf, da_snowc, rain_threshold)
        da_lists["rain_on_snow_freezing"].append(rosf_tmp.compute())
        da_lists["rain_on_snow_thawing"].append(rost_tmp.compute())
        
    # # winter warming events 
    if variables["winter_warming_events"]:
        wwe_tmp = indices.winter_warming_events(da_t2mean_winter, da_snowc)
        da_lists["winter_warming_events"].append(wwe_tmp.compute())

    # # winter warming events 
    if variables["winter_warming_intensity"]:
        wwi_tmp = indices.winter_warming_intensity(da_t2mean_winter, da_snowc)
        da_lists["winter_warming_intensity"].append(wwi_tmp.compute())

    # # Heatwave magnitude index
    if variables["heatwave_magnitude"]:
        hwmi_tmp = indices.heatwave_magnitude_index(da_t2max, T90p, p75max, p25max)
        da_lists["heatwave_magnitude"].append(hwmi_tmp.compute())

    # # Vapour pressure deficit index   
    if variables["vpd_magnitude_index"]:    
        vpdi_tmp = indices.vpd_magnitude_index(da_t2mean_summer, da_d2mean, T90p_vpd, p75vpd, p25vpd)
        da_lists["vpd_magnitude_index"].append(vpdi_tmp.compute())
  
    # # Summer warmth index
    if variables["summer_warmth_index"]:
        swi_tmp = indices.summer_warmth_index(da_t2mean_summer)
        da_lists["summer_warmth_index"].append(swi_tmp.compute())
        
    # # Gale wind events
    if variables["high_wind_events"]:
        hwe_tmp = indices.high_wind_events(da_u10max, da_v10max)
        da_lists["high_wind_events"].append(hwe_tmp.compute())
    
    # # Annual mean temperature
    if variables["annual_mean_temperature"]:
        tavg_tmp = indices.annual_mean_temperature(da_t2mean_summer)
        da_lists["annual_mean_temperature"].append(tavg_tmp.compute())

    # # Annual precipitation
    if variables["annual_precipitation"]:
        tpa_tmp = indices.annual_precipitation(da_tp_summer)
        da_lists["annual_precipitation"].append(tpa_tmp.compute())

    # # Annual snowfall
    if variables["annual_snowfall"]:
        sfa_tmp = indices.annual_snowfall(da_sf)
        da_lists["annual_snowfall"].append(sfa_tmp.compute())
        
    # # Average wind speed
    if variables["average_wind_speed"]:
        aws_tmp = indices.average_wind_speed(da_u10mean, da_v10mean)
        da_lists["average_wind_speed"].append(aws_tmp.compute())
        
    
    # print memory usage
    print('Memory usage:')
    print(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, 'MB')

print('All incides computed!')

# concatenate dataarrays and save data into netcdf files
for var in da_lists:
    ds_out[var] = xr.concat(da_lists[var], dim='time').compute().to_dataset(name=shortnames[var])
    
    # add time attribute
    ds_out[var].time.attrs['long_name'] = "time"

    # add global attributes
    ds_out[var].attrs['Conventions'] = 'CF-1.7'
    ds_out[var].attrs['title'] = 'ARCLIM Bioclimatic indices'
    ds_out[var].attrs['Institution'] = 'Finnish Meteorological Institute'
    ds_out[var].attrs['source'] = 'ERA5-Land'
    ds_out[var].attrs['history'] = datetime.utcnow().strftime(format='%Y-%m-%d %H:%M:%S') + ' Python'
    
    # define outfile
    outfile_nc = outpath  + '/nc/arclim_' + shortnames[var] + '.nc'
    outfile_tiff = outpath  + '/tif/arclim_' + shortnames[var] + '.tif'
    
    # save the data as a netcdf file
    ds_out[var].to_netcdf(outfile_nc, format='NETCDF4', encoding={'time': {'dtype': 'i4'}})
    
    # convert nc to geotiff and save
    geotiff = io_utils.da_to_geotiff(ds_out[var][shortnames[var]])
    geotiff.to_raster(outfile_tiff)

    
    


print('Done!')
