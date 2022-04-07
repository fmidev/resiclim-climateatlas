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
import ghcn_routines as ghcn
from scipy import stats
import math

def rounddown(x):
    return int(math.floor(x / 100.0)) * 100

def roundup(x):
    return int(math.ceil(x / 100.0)) * 100



# list of stations and their names
list_of_stations = ghcn.ghcn_stations()

abcs = ['a)','b)','c)','d)','e)','f)', 'g)', 'h)','i)']


# read ARCLIM dataset
arclim_ds = xr.open_dataset('/Users/rantanem/Downloads/arclim_FDD.nc')
fdd_da = arclim_ds['FDD']


# read station locations
locs = ghcn.read_station_locations()

# read the station data
fdd = pd.read_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_fdd.csv',
                  index_col=0)
year1=fdd.index[0]
year2=fdd.index[-1]


# read station location coordinates from GHCN server
station_locs = pd.DataFrame(index=fdd.columns, columns=['lat','lon'])
for station in fdd.columns:
    station_locs.loc[station].lat = locs.loc[station].lat
    station_locs.loc[station].lon = locs.loc[station].lon



fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), constrained_layout=False, dpi=200,)

axlist = axarr.flatten()

plt.subplots_adjust(wspace=0.4, hspace=0.4)

for i, station in enumerate(fdd.columns):
    
    lat = station_locs.iloc[i].lat
    lon = station_locs.iloc[i].lon
    
    
    station_fdd = fdd[station_locs.iloc[i].name]
    idx = np.isfinite(station_fdd.values)
    station_fdd = station_fdd[idx]

    arclim_fdd = fdd_da.sel(latitude=lat, longitude=lon, method='nearest').sel(time=slice(year1, year2))
    
    arclim_fdd = arclim_fdd[idx]
    
    x = station_fdd
    y = arclim_fdd.values
    
    axlist[i].scatter(x, y, color='lightseagreen')
    axlist[i].plot([0, 1], [0, 1], transform=axlist[i].transAxes, color='lightseagreen')

    result = stats.linregress(station_fdd,arclim_fdd.values)
    r = np.round(result.rvalue, 2)
    bias = np.round(np.mean(arclim_fdd.values - station_fdd))
    
    axlist[i].annotate(abcs[i], (0.04, 0.91), xycoords='axes fraction', ha='left', 
                fontweight='bold', fontsize=18)
    axlist[i].annotate('R: '+str(r), (0.93, 0.03), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    axlist[i].annotate('Bias: '+str(bias), (0.93, 0.1), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    
    
    xmin = rounddown(np.min([arclim_fdd.min().values-100, station_fdd.min()-100]))
    xmax = roundup(np.max([arclim_fdd.max().values+100, station_fdd.max()+100]))
    # axlist[i].set_xticks(np.arange(xmin, xmax+600, 600))
    # axlist[i].set_yticks(np.arange(xmin, xmax+600, 600))
    axlist[i].set_xlim(xmin, xmax)
    axlist[i].set_ylim(xmin, xmax)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    axlist[i].set_ylabel('ARCLIM FDD [°C days]', fontsize=14)
    axlist[i].set_xlabel('Station FDD [°C days]', fontsize=14)
    
    axlist[i].set_title(list_of_stations[station], fontsize=16)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_scatterplot_fdd.png'
   
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')



fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), constrained_layout=False, dpi=200,)

axlist = axarr.flatten()

for i, station in enumerate(fdd.columns):
    
    lat = station_locs.iloc[i].lat
    lon = station_locs.iloc[i].lon
    
    
    station_fdd = fdd[station_locs.iloc[i].name]
    idx = np.isfinite(station_fdd.values)
    # station_gdd = station_gdd[idx]
   
    arclim_fdd = fdd_da.sel(latitude=lat, longitude=lon, method='nearest').sel(time=slice(year1, year2))
    # arclim_gdd = arclim_gdd[idx]
    
    
    axlist[i].plot(arclim_fdd.time, arclim_fdd.values, '-o', color='lightseagreen', label='ARCLIM')
    axlist[i].plot(station_fdd.index, station_fdd, '-o', color='darkorange', label='Station')

    result = stats.linregress(station_fdd[idx],arclim_fdd[idx].values)
    rsquared = np.round(result.rvalue**2, 2)
    
    axlist[i].annotate('R²: '+str(rsquared), (0.04, 0.92), xycoords='axes fraction', ha='left', 
                fontstyle='italic', fontsize=14)
    
    axlist[i].set_xlim(1959, 2022)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    
    axlist[i].set_title(list_of_stations[station], fontsize=14)    
    
axlist[0].legend(bbox_to_anchor=(2.24,1.3), ncol=2, fontsize=18)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_timeseries_fdd.png'  
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')
   
