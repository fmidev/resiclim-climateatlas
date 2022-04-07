#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 08:34:42 2022

These routines read data from GHCN-D database. The stations used for ARCLIM 
validation are defined here.

@author: rantanem
"""
import pandas as pd
import numpy as np

def ghcn_stations():
    
    # list of stations and their names
    list_of_stations = {'SV000001008':'Svalbard, Norway',
                        'FI000007501':'Sodankylä, Finland', 
                        'RSM00027612':'Moscow, Russia',
                        'RSM00021802':'Saskylah, Russia',
                        'RSM00024266':'Verhojansk, Russia',
                        'USW00025503':'King Salmon, US',
                        'CA002300500':'Baker Lake, Canada',
                        'IC000004030':'Reykjavik, Iceland',
                        'CA002401200':'Eureka, Canada',}
    
    return list_of_stations

def ghcn_station_offsets():
    
    # these offsets are used in plotting the station names to the map
    
    offsets = {'SV000001008':(28,-10),
               'FI000007501':(28,-10),
               'RSM00021802':(24,9),
               'RSM00024266':(24,9),
               'RSM00027612':(24,-11),
               'USW00025503':(38,9),
               'CA002300500':(-26,9),
               'IC000004030':(-24,-11),
               'CA002401200':(24,9)}
    
    return offsets

def read_station_locations():
    
    ## read all station locations
    locs = pd.read_csv('https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt', 
                       delim_whitespace=True, 
                       names=['station','lat','lon','D','E','F','G','H','I','J','K','L','M'])
    locs.index = locs.station
    
    stations = ghcn_stations()
    
    # select only the chosen stations
    station_locs = pd.DataFrame(index=stations.keys(), columns=['lat','lon'])
    for station in stations.keys():
        station_locs.loc[station].lat = locs.loc[station].lat
        station_locs.loc[station].lon = locs.loc[station].lon
    
    return station_locs

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
