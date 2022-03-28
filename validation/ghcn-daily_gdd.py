#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 11:18:02 2022

@author: rantanem
"""
import pandas as pd
import numpy as np
from routines import update_station_data


def get_ghcn_daily_var(var, station, years):
    
    url = 'https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/'+station+'.csv.gz'
    
    colnames=['staid', 'date', 'variable', 'value', 'flag1', 'flag2', 'flag3','flag4']
    
    ghcn_daily = pd.read_csv(url,low_memory=False,
                         names=colnames, header=None)
    
    f = ghcn_daily[ghcn_daily.variable==var]
    f.index = pd.to_datetime(f.date, format='%Y%m%d')
    f = f[['value']]
    f['value'] = f['value']/10
    f.rename(columns={'value':var}, inplace=True)
    
    cond = np.isin(f.index.year, years)
    
    
    return f[cond]


list_of_stations = ['FI000007501', 'RSM00023804', 'RSM00024136',
                    'RSM00025551', 'IC000004030', 'CA002300500', 'USW00025503',]



basevalue = 5

years = np.arange(1960,2022)
dates = pd.date_range(str(years[0])+'-01-01', str(years[-1])+'-12-31')

df_daily_data = pd.DataFrame(index=dates, columns=list_of_stations)
df_gdd = pd.DataFrame(index=years, columns=list_of_stations)



# get the data; loop over stations
for i, station in enumerate(list_of_stations):
    
    print(station)
    
    if station[:2]=='FI':
        dataset, sname = update_station_data(station='sodankylä', forecast=False)
        cond = np.isin(dataset.index.year, years)
        f = dataset['Average temperature'][cond]
    elif (station[:2]=='US') | (station[:2]=='CA'):
        f1 = get_ghcn_daily_var('TMAX', station, years)
        f2 = get_ghcn_daily_var('TMIN', station, years)
        f = (f1.TMAX+f2.TMIN)/2
        
    else: 
        f = get_ghcn_daily_var('TAVG', station, years)
    
    df_daily_data[station] = f.reindex(dates)
    print('Number of missing values:',np.sum(f.reindex(dates).isna()))
    
# calculate growing degree days
for i, station in enumerate(list_of_stations):
    
    station_data = df_daily_data[station]
    
    #loop through spring
    for y in years:
        cond = (station_data.index.year == y)  & (station_data.index.month  > 0)
        data = station_data[cond]
        # data.index = pd.to_datetime(station_data[cond][['Year','Month','Day']])
        temp = data - basevalue
        N = np.sum(temp[str(y)+'-04-01':str(y)+'-09-30'].isna())
        # print(y, N)
        
        if N > 2: # check if there are too much missing values
            gsbeg_idx = pd.NaT
            cs = np.nan
        elif len(temp) == 0: # check if the data is missing
            gsbeg_idx = pd.NaT
            cs = np.nan
        else:
            gsbeg_idx = temp[str(y)+'-01-01':str(y)+'-06-30'].cumsum(skipna=True).idxmin() + pd.Timedelta(1, unit = 'D')
            gsend_idx = temp[str(y)+'-07-01':str(y)+'-12-31'].cumsum(skipna=True).idxmax()

            gdd = temp[gsbeg_idx:gsend_idx]
        
            gdd[gdd<0] = 0
        
            df_gdd[station][y] = gdd.sum(skipna=True)
    
    
df_gdd.to_csv('/Users/rantanem/Documents/python/data/ghcn_daily_gdd.csv',)

df_gdd.plot()