#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 14:11:12 2022

This script compares observed and ERA5-Land summer average temperatures. 
The observed temperatures are read from GHCN-Daily database.
The results are plotted as a 9-panel scatterplot 

@author: rantanem
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
import math
import ghcn_routines as ghcn

def rounddown(x):
    return int(math.floor(x / 10.0)) * 10

def roundup(x):
    return int(math.ceil(x / 10.0)) * 10


# list of stations and their names
list_of_stations = ghcn.ghcn_stations()

abcs = ['a)','b)','c)','d)','e)','f)', 'g)', 'h)','i)']

# read ARCLIM dataset
arclim_ds = xr.open_dataset('/Users/rantanem/Downloads/resiclim_annual.nc')
gdd_da = arclim_ds['Tavg_summer']-273.15


# read station locations
station_locs = ghcn.read_station_locations()

# read the station data
gdd = pd.read_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/stations_daily_tavg_summer.csv',
                  index_col=0)
year1=gdd.index[0]
year2=gdd.index[-1]



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
    bias = np.round(np.mean(arclim_gdd.values - station_gdd),1)
    
    axlist[i].annotate(abcs[i], (0.04, 0.92), xycoords='axes fraction', ha='left', 
                fontweight='bold', fontsize=18)
    axlist[i].annotate('R: '+str(r), (0.95, 0.03), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    axlist[i].annotate('Bias: '+str(bias), (0.95, 0.1), xycoords='axes fraction', ha='right', 
                fontstyle='italic', fontsize=14)
    
    
    xmin = np.min([arclim_gdd.min().values-1, station_gdd.min()-1])
    xmax = np.max([arclim_gdd.max().values+1, station_gdd.max()+1])
    # axlist[i].set_xticks(np.arange(xmin, xmax+200, 200))
    # axlist[i].set_yticks(np.arange(xmin, xmax+200, 200))
    axlist[i].set_xlim(xmin, xmax)
    axlist[i].set_ylim(xmin, xmax)
    axlist[i].tick_params(axis='both', which='major', labelsize=14)
    axlist[i].tick_params(axis='both', which='minor', labelsize=14)
    axlist[i].set_ylabel('ERA5-Land [°C]', fontsize=14)
    axlist[i].set_xlabel('Station [°C]', fontsize=14)
    
    axlist[i].set_title(list_of_stations[station_locs.iloc[i].name], fontsize=16)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_scatterplot_tavg_summer.png'
   
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')