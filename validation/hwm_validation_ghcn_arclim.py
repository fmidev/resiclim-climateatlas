#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:33:54 2022

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
    return int(math.floor(x / 10.0)) * 10

def roundup(x):
    return int(math.ceil(x / 10.0)) * 10



# list of stations and their names
list_of_stations = ghcn.ghcn_stations()

abcs = ['a)','b)','c)','d)','e)','f)', 'g)', 'h)','i)']

# read ARCLIM dataset
arclim_ds = xr.open_dataset('/Users/rantanem/Downloads/arclim_HWM.nc')
hwm_da = arclim_ds['HWM']# - 273.15


# read station locations
locs = ghcn.read_station_locations()

# read the station data
hwm = pd.read_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/stations_daily_hwm.csv',
                  index_col=0)
year1=hwm.index[0]
year2=hwm.index[-1]


# read station location coordinates from GHCN server
station_locs = pd.DataFrame(index=hwm.columns, columns=['lat','lon'])
for station in hwm.columns:
    station_locs.loc[station].lat = locs.loc[station].lat
    station_locs.loc[station].lon = locs.loc[station].lon



fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), constrained_layout=False, dpi=200,)

axlist = axarr.flatten()

plt.subplots_adjust(wspace=0.4, hspace=0.4)

for i, station in enumerate(hwm.columns):
    
    lat = station_locs.iloc[i].lat
    lon = station_locs.iloc[i].lon
    
    
    station_hwm = hwm[station_locs.iloc[i].name]
    idx = np.isfinite(station_hwm.values)
    station_hwm = station_hwm[idx]

    arclim_hwm = hwm_da.sel(latitude=lat, longitude=lon, method='nearest').sel(time=slice(year1, year2))
    
    arclim_hwm = arclim_hwm[idx]
    
    # top_years = station_hwm.sort_values().index[-1:] 
    # idx = np.isin(station_hwm.index, top_years, invert=True)
    x = station_hwm
    y = arclim_hwm.values
    
    axlist[i].scatter(x, y, color='lightseagreen')
    axlist[i].plot([0, 1], [0, 1], transform=axlist[i].transAxes, color='lightseagreen')


    result = stats.linregress(x,y)
    r = np.round(result.rvalue, 2)
    
    result = stats.spearmanr(x, y)
    r_spear = np.round(result.correlation, 2)
    
    axlist[i].annotate(abcs[i], (0.04, 0.91), xycoords='axes fraction', ha='left', 
                fontweight='bold', fontsize=18)
    axlist[i].annotate('R: '+str(r_spear), (0.93, 0.03), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    
    bias = np.round(arclim_hwm.mean().values - station_hwm.mean(), 2)
    axlist[i].annotate('Bias: '+str(bias), (0.93, 0.1), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    
    

    xmax = np.ceil(np.max([arclim_hwm.max().values+1, station_hwm.max()+1]))
    # axlist[i].set_xticks(np.arange(xmin, xmax+600, 600))
    # axlist[i].set_yticks(np.arange(xmin, xmax+600, 600))
    axlist[i].set_xlim(-1, xmax)
    axlist[i].set_ylim(-1, xmax)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    axlist[i].set_ylabel('ARCLIM HWM', fontsize=14)
    axlist[i].set_xlabel('Station HWM', fontsize=14)
    
    axlist[i].set_title(list_of_stations[station], fontsize=16)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_scatterplot_hwm.png'
   
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')



fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), constrained_layout=False, dpi=200,)

axlist = axarr.flatten()

for i, station in enumerate(hwm.columns):
    
    lat = station_locs.iloc[i].lat
    lon = station_locs.iloc[i].lon
    
    
    station_hwm = hwm[station_locs.iloc[i].name]
    idx = np.isfinite(station_hwm.values)
    # station_gdd = station_gdd[idx]
   
    arclim_hwm = hwm_da.sel(latitude=lat, longitude=lon, method='nearest').sel(time=slice(year1, year2))
    # arclim_gdd = arclim_gdd[idx]
    
    
    axlist[i].plot(arclim_hwm.time, arclim_hwm.values, '-o', color='lightseagreen', label='ARCLIM')
    axlist[i].plot(station_hwm.index, station_hwm, '-o', color='darkorange', label='Station')

    result = stats.linregress(station_hwm[idx],arclim_hwm[idx].values)
    rsquared = np.round(result.rvalue**2, 2)
    
    axlist[i].annotate('R²: '+str(rsquared), (0.04, 0.92), xycoords='axes fraction', ha='left', 
                fontstyle='italic', fontsize=14)
    
    # axlist[i].set_xlim(1980, 2011)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    
    axlist[i].set_title(list_of_stations[station], fontsize=14)    
    
axlist[0].legend(bbox_to_anchor=(2.24,1.3), ncol=2, fontsize=18)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_timeseries_hwm.png'  
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')
   
