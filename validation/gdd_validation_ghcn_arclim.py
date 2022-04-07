#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 14:11:12 2022

@author: rantanem
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats, signal

import math
from ghcn_routines import ghcn_stations

def rounddown(x):
    return int(math.floor(x / 100.0)) * 100

def roundup(x):
    return int(math.ceil(x / 100.0)) * 100

def read_station_locations():
    
    ## read station location
    locs = pd.read_csv('https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt', 
                       delim_whitespace=True, 
                       names=['station','lat','lon','D','E','F','G','H','I','J','K','L','M'])
    locs.index = locs.station
    
    return locs[['lat','lon']]


# list of stations and their names
list_of_stations = ghcn_stations()

abcs = ['a)','b)','c)','d)','e)','f)', 'g)', 'h)','i)']

# read ARCLIM dataset
arclim_ds = xr.open_dataset('/Users/rantanem/Downloads/resiclim_GDD.nc')
gdd_da = arclim_ds['GDD']


# read station locations
locs = read_station_locations()

# read the station data
gdd = pd.read_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_gdd.csv',
                  index_col=0)
year1=gdd.index[0]
year2=gdd.index[-1]


# read station location coordinates from GHCN server
station_locs = pd.DataFrame(index=gdd.columns, columns=['lat','lon'])
for station in gdd.columns:
    station_locs.loc[station].lat = locs.loc[station].lat
    station_locs.loc[station].lon = locs.loc[station].lon



fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), constrained_layout=False, dpi=200,)

axlist = axarr.flatten()

plt.subplots_adjust(wspace=0.4, hspace=0.4)

for i, ax in enumerate(axlist):
    
    lat = station_locs.iloc[i].lat
    lon = station_locs.iloc[i].lon
    
    
    station_gdd = gdd[station_locs.iloc[i].name]
    idx = np.isfinite(station_gdd.values)
    station_gdd = station_gdd[idx]

    arclim_gdd = gdd_da.sel(latitude=lat, longitude=lon, method='nearest').sel(time=slice(year1, year2))
    
    arclim_gdd = arclim_gdd[idx]
    
    axlist[i].scatter(station_gdd, arclim_gdd, color='lightseagreen')
    axlist[i].plot([0, 1], [0, 1], transform=axlist[i].transAxes, color='lightseagreen')

    x = station_gdd
    y = arclim_gdd.values
    
    result = stats.linregress(x,y)
    r = np.round(result.rvalue, 2)
    bias = np.round(np.mean(arclim_gdd.values - station_gdd), 1)
    
    axlist[i].annotate(abcs[i], (0.04, 0.92), xycoords='axes fraction', ha='left', 
                fontweight='bold', fontsize=18)
    axlist[i].annotate('R: '+str(r), (0.93, 0.03), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    axlist[i].annotate('Bias: '+str(bias), (0.93, 0.1), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    
    
    xmin = rounddown(np.min([arclim_gdd.min().values-10, station_gdd.min()-10]))
    xmax = roundup(np.max([arclim_gdd.max().values+10, station_gdd.max()+10]))
    # axlist[i].set_xticks(np.arange(xmin, xmax+200, 200))
    # axlist[i].set_yticks(np.arange(xmin, xmax+200, 200))
    axlist[i].set_xlim(xmin, xmax)
    axlist[i].set_ylim(xmin, xmax)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    axlist[i].set_ylabel('ARCLIM GDD [°C days]', fontsize=14)
    axlist[i].set_xlabel('Station GDD [°C days]', fontsize=14)
    
    axlist[i].set_title(list_of_stations[station_locs.iloc[i].name], fontsize=16)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_scatterplot_gdd.png'
   
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')



fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), constrained_layout=False, dpi=200,)

axlist = axarr.flatten()

for i, ax in enumerate(axlist):
    
    lat = station_locs.iloc[i].lat
    lon = station_locs.iloc[i].lon
    
    
    station_gdd = gdd[station_locs.iloc[i].name]
    idx = np.isfinite(station_gdd.values)
    # station_gdd = station_gdd[idx]
   
    arclim_gdd = gdd_da.sel(latitude=lat, longitude=lon, method='nearest').sel(time=slice(year1, year2))
    # arclim_gdd = arclim_gdd[idx]
    
    
    axlist[i].plot(arclim_gdd.time, arclim_gdd.values, '-o', color='lightseagreen', label='ARCLIM')
    axlist[i].plot(station_gdd.index, station_gdd, '-o', color='darkorange', label='Station')

    # result = stats.linregress(station_gdd[idx],arclim_gdd[idx].values)
    # rsquared = np.round(result.rvalue**1, 2)
    
    x = station_gdd[idx].index
    y = station_gdd[idx]
    
    res = stats.theilslopes(y, x, 0.90)
    ax.plot(x, res[1] + res[0] * x, linestyle='-', color='darkorange')
    trend = int(np.round(res[0]*10,0))
    
    axlist[i].annotate('Trend: '+str(trend) + ' °C days dec⁻¹', (0.15, 0.93), xycoords='axes fraction', ha='left', 
                fontstyle='italic', fontsize=14, color='darkorange')
    
    x = station_gdd[idx].index
    y = arclim_gdd[idx]
    
    res = stats.theilslopes(y, x, 0.90)
    ax.plot(x, res[1] + res[0] * x, linestyle='-', color='lightseagreen')
    trend = int(np.round(res[0]*10,0))
    
    axlist[i].annotate('Trend: '+str(trend) + ' °C days dec⁻¹', (0.15, 0.87), xycoords='axes fraction', ha='left', 
                fontstyle='italic', fontsize=14, color='lightseagreen')
    
    axlist[i].annotate(abcs[i], (0.04, 0.92), xycoords='axes fraction', ha='left', 
                fontweight='bold', fontsize=18)
    xmax = roundup(np.max([arclim_gdd.max().values+100, station_gdd.max()+100]))
    xmin = rounddown(np.min([arclim_gdd.min().values-10, station_gdd.min()-10]))
    axlist[i].set_ylim(xmin, xmax)
    axlist[i].set_xlim(1959, 2022)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    
    axlist[i].set_title(list_of_stations[station_locs.iloc[i].name], fontsize=16)
    
axlist[0].legend(bbox_to_anchor=(2.24,1.3), ncol=2, fontsize=18)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_timeseries_gdd.png'  
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')
   
