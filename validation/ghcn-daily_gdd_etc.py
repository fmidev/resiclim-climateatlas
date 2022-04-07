#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 11:18:02 2022

This script reads daily temperature values from NOAA GHCN-Daily archive,
and calculates growing season degree days, freezing degree days and 
summer mean temperatures. FMI Sodankylä station data is read from FMI archives. 

For Russian and Icelandic stations, TAVG (daily average temperature) is read.
For North American stations, TMAX and TMIN (daily max and min temperatures) 
are read, and TAVG is calculated as the average of these two.

@author: rantanem
"""
import pandas as pd
import numpy as np
from fmi_routines import update_station_data
import ghcn_routines as ghcn


# list of stations and their names
list_of_stations = ghcn.ghcn_stations()

# base value for thermal growing degree days in celcius
basevalue = 5

# years and dates for which the GDD is calculated
years = np.arange(1960,2022)
dates = pd.date_range(str(years[0])+'-01-01', str(years[-1])+'-12-31')

# allocate dataframes
df_daily_data = pd.DataFrame(index=dates, columns=list_of_stations)
df_gdd = pd.DataFrame(index=years, columns=list_of_stations)
df_fdd = pd.DataFrame(index=years, columns=list_of_stations)
df_tavg_summer = pd.DataFrame(index=years, columns=list_of_stations)



# get the data; loop over stations
for i, station in enumerate(list_of_stations):
    
    print(list_of_stations[station])
    
    # Finnish data is read from FMI
    if station[:2]=='FI':
        dataset = update_station_data(station='sodankyla')
        cond = np.isin(dataset.index.year, years)
        f = dataset['Average temperature'][cond]
    # US and Canadian TAVG is calculated from TX and TN
    elif (station[:2]=='US') | (station[:2]=='CA'):
        f1 = ghcn.get_ghcn_daily_var('TMAX', station, years)
        f2 = ghcn.get_ghcn_daily_var('TMIN', station, years)
        f = (f1.TMAX+f2.TMIN)/2
    # for other stations, TAVG is read    
    else: 
        f = ghcn.get_ghcn_daily_var('TAVG', station, years)
    
    # allocate data to the dataframe
    df_daily_data[station] = f.reindex(dates)
    
    print('Number of missing values:',np.sum(f.reindex(dates).isna().values),'\n')
  
    

# calculate the growing degree days
for i, station in enumerate(list_of_stations):
    
    station_data = df_daily_data[station]
    
    #loop through the years
    for y in years:
        cond = (station_data.index.year == y)  & (station_data.index.month  > 0)
        data = station_data[cond]

        # subtract the basevalue
        temp = data - basevalue
        # number of missing days in Apr-Sep period
        N = np.sum(temp[str(y)+'-04-01':str(y)+'-09-30'].isna())

        
        if N > 1: # check if there are too much missing values
            gsbeg_idx = pd.NaT
            cs = np.nan
        elif len(temp) == 0: # check if the data is missing
            gsbeg_idx = pd.NaT
            cs = np.nan
        else:
            gsbeg_idx = temp[str(y)+'-01-01':str(y)+'-06-30'].cumsum(skipna=True).idxmin() + pd.Timedelta(1, unit = 'D')
            gsend_idx = temp[str(y)+'-07-01':str(y)+'-12-31'].cumsum(skipna=True).idxmax()

            gdd = temp[gsbeg_idx:gsend_idx]
        
            # negative degrees are marked as 0
            gdd[gdd<0] = 0
        
            df_gdd[station][y] = gdd.sum(skipna=True)
    

# save the GDD values
df_gdd.to_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_gdd.csv',)

df_gdd.plot()

# calculate the freezing degree days
for i, station in enumerate(list_of_stations):
    
    station_data = df_daily_data[station]
    
    #loop through the years
    for y in years[1:]:
        cond1 = (station_data.index.year == y-1)  & (station_data.index.month >7)
        cond2 = (station_data.index.year == y)  & (station_data.index.month <8)
        temp = station_data[cond1+cond2]

        # number of missing days in Apr-Sep period
        N = np.sum(temp[str(y-1)+'-10-01':str(y)+'-04-30'].isna())

        
        if N > 1: # check if there are too much missing values
            beg_idx = pd.NaT
            cs = np.nan
        elif len(temp) == 0: # check if the data is missing
            beg_idx = pd.NaT
            cs = np.nan
        else:
            beg_idx = temp[str(y-1)+'-08-01':str(y)+'-01-31'].cumsum(skipna=True).idxmax() + pd.Timedelta(1, unit = 'D')
            end_idx = temp[str(y)+'-02-01':str(y)+'-07-31'].cumsum(skipna=True).idxmin()

            fdd = temp[beg_idx:end_idx]
        
            # negative degrees are marked as 0
            fdd[fdd>0] = 0
            
            # multiply by -1 to get positive values
            fdd *= -1
        
            df_fdd[station][y] = fdd.sum(skipna=True)
    

# save the FDD values
df_fdd.to_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_fdd.csv',)

df_fdd.plot()

# calculate the average summer temperature
for i, station in enumerate(list_of_stations):
    
    station_data = df_daily_data[station]
    
    #loop through the years
    for y in years:
        cond = (station_data.index.year == y)  & (np.isin(station_data.index.month, [6,7,8]))
        temp = station_data[cond]
        
        # number of missing days in Apr-Sep period
        N = np.sum(temp.isna())
        
        if N > 1: # check if there are too much missing values
            value=np.nan
        else:
            value=temp.mean()
        
            df_tavg_summer[station][y] = value

df_tavg_summer.to_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_tavg_summer.csv',)
