#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 11:18:02 2022

This script reads daily temperature values from NOAA GHCN-Daily archive,
and calculates heatwave magnitude index. FMI Sodankylä station data is 
read from FMI database. 

@author: rantanem
"""
import pandas as pd
import numpy as np
from fmi_routines import update_station_data
import ghcn_routines as ghcn
from scipy.ndimage.morphology import binary_dilation
from scipy import ndimage

# take only summer months
def is_jja(month):
    return (month >= 6) & (month <= 8)

# list of stations and their names
list_of_stations = ghcn.ghcn_stations()

# read station location coordinates from GHCN server
station_locs = ghcn.read_station_locations()

# years and dates for which the HWM is calculated
years = np.arange(1960,2022)
dates = pd.date_range(str(years[0])+'-01-01', str(years[-1])+'-12-31')

# allocate empty dataframes
df_daily_data = pd.DataFrame(index=dates, columns=list_of_stations)
df_hwmi = pd.DataFrame(index=years, columns=list_of_stations)
df_tmax = pd.DataFrame(index=years, columns=list_of_stations)



# get the data; loop over the stations
for i, station in enumerate(list_of_stations):
    
    print(list_of_stations[station])
    
    # Finnish data is read from FMI
    if station[:2]=='FI':
        dataset = update_station_data(station='sodankyla')
        cond = np.isin(dataset.index.year, years)
        f = dataset['Maximum temperature'][cond]
    # for other stations, read TX from GHCN-Daily
    else:
        f = ghcn.get_ghcn_daily_var('TMAX', station, years)

    # allocate data to the dataframe
    df_daily_data[station] = f.reindex(dates)
    
    # print the number of missing days
    print('Number of missing values:',np.sum(f.reindex(dates).isna().values),'\n')
  

# Width of the threshold selection window (days)
struct = np.ones(31) 
df_p90 = pd.DataFrame(index=np.unique(df_daily_data.index.dayofyear), columns=list_of_stations)
df_25_75 = pd.DataFrame(index=[25, 75], columns=list_of_stations)

# climatology years for the threshold 
years_clim = np.arange(1981, 2011)


# calculate the threshold for heat wave magnitude index
# (the 90th percentile of daily maximum temperature)
for i, station in enumerate(list_of_stations):
    
    station_data_all_years = df_daily_data[station]
    
    # select only the 1981-2010 years
    cond = np.isin(station_data_all_years.index.year, years_clim)
    station_data = station_data_all_years[cond]
    
    doy_values = np.unique(station_data.index.dayofyear)
    
    # Loop over each day of year
    for day in doy_values:
        
        dayofyear = station_data.index.dayofyear == day
        
        selection = binary_dilation(dayofyear, structure=struct)
        
        temp = station_data[selection]
        
        df_p90[station][day] = np.nanpercentile(temp, 90)
 

# calculate the 25th and 75th percentiles of annual maxima
for i, station in enumerate(list_of_stations):
    
    station_data = df_daily_data[station]
    years_clim = np.arange(1981, 2011)
    cond = np.isin(station_data.index.year, years_clim)
    station_data = station_data[cond]
    
    maxvalues = station_data.groupby(station_data.index.year).max()
    
    p75_max = maxvalues.quantile(0.75)
    p25_max = maxvalues.quantile(0.25)
    
    df_25_75[station][25] = p25_max
    df_25_75[station][75] = p75_max
    



# generate the structure to label each heatwave event
struct = np.ones(shape=(3,))

    
# calculate the heat wave magnitude index
for i, station in enumerate(list_of_stations):    
    
    station_data = df_daily_data[station]
    heatwave_threshold = df_p90[station]
    #loop through the years
    for y in years:
        cond = (station_data.index.year == y)  & (station_data.index.month  > 0)
        temp = station_data[cond]
        N = np.sum(temp[str(y)+'-06-01':str(y)+'-08-31'].isna())
        
        newcoords = pd.to_datetime(y * 1000 + df_p90.index, format='%Y%j')   
        
        heatwave_threshold.index = newcoords
        
        # identify heatwave days
        heatwaves = temp > heatwave_threshold[temp.index]
        
        # label each heatwave event
        labels, nb = ndimage.label(heatwaves, structure=struct)
        
        # calculate the length of each heatwave
        heatwave_lengths = np.array(ndimage.sum(heatwaves, labels, np.arange(labels.max()+1)))
        
        # mask heatwaves which are shorther than three days
        mask = heatwave_lengths > 2
        heatwave_events = mask[labels.ravel()].reshape(labels.shape)

        # select only JJA period
        heatwave_temps =  temp[heatwave_events][is_jja(temp[heatwave_events].index.month)]
        heatwave_temps =  temp.where(heatwave_events)[is_jja(temp.where(heatwave_events).index.month)]
        
        # calculate the heatwave magnitude (based on Dobricic et al. 2020)
        Md = (heatwave_temps - df_25_75[station][25]) / (df_25_75[station][75] - df_25_75[station][25])

        # calculate sum of the daily magnitudes of Md from the consecutive days composing a heat wave
        # negative Md indices are considered zero
        heatwavevalues = (Md.where((Md>0)|(Md.isnull()), 0)).values
        cums = np.nancumsum(heatwavevalues, axis=0)
        weights_by_duration_array = cums - np.maximum.accumulate(cums * (np.isnan(heatwavevalues)), axis=0)

        hwi = np.max(weights_by_duration_array)

        if N == 0:
            df_hwmi[station][y] = hwi
            df_tmax[station][y] = temp.max()

# save the HWM values
df_hwmi.to_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_hwm.csv',
               index_label='Year', na_rep='NaN')
df_tmax.to_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_tmax.csv',
               index_label='Year', na_rep='NaN')
  