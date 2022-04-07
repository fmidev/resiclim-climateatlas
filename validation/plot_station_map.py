#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 12:49:08 2022

@author: rantanem
"""
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import ghcn_routines

def read_station_locations():
    
    ## read station location
    locs = pd.read_csv('https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt', 
                       delim_whitespace=True, 
                       names=['station','lat','lon','D','E','F','G','H','I','J','K','L','M'])
    locs.index = locs.station
    
    return locs[['lat','lon']]

def plotMap():


    #Set the projection information, depending on the area you want to plot
    proj = ccrs.NorthPolarStereo()
    #Create a figure with an axes object on which we will plot. Pass the projection to that axes.
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(11, 6), constrained_layout=False, dpi=200,
                          subplot_kw={'projection': proj})
    # plot grid lines if needed
    # gl = ax.gridlines(linewidth=2, color='k', alpha=0.5, linestyle='--')
    # gl.n_steps = 100
    # # define which grid lines you want to plot
    # gl.ylocator = mticker.FixedLocator([60, 75])
    # gl.xlocator = mticker.FixedLocator([60, 120])
    plot_background(ax)
   
    return fig, ax

def plot_background(ax):
    import cartopy.feature as cfeature
    import matplotlib.path as mpath  

    ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
    # ax.set_extent([-20, 45, 33,74], crs=ccrs.PlateCarree())
    # ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor=cfeature.COLORS['land'])
    # ax.add_feature(cfeature.LAKES.with_scale('50m'),facecolor=cfeature.COLORS['land'],zorder=1,edgecolor='k',linewidth=1.5) 
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), zorder=5, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), zorder=5,edgecolor='k',linewidth=0.5)
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)
    
    # gl = ax.gridlines(linewidth=2, color='k', alpha=0.5, linestyle='--')
    # gl.n_steps = 100
    # gl.ylocator = mticker.FixedLocator([66.5])
    # gl.xlocator = mticker.FixedLocator([85, 160, 160])
    return ax


# read station locations
locs = read_station_locations()

# list of stations and their names
list_of_stations = ghcn_routines.ghcn_stations()

# list of stations and their names
offsets = ghcn_routines.ghcn_station_offsets()

station_locs = pd.DataFrame(index=list(list_of_stations.keys()), columns=['lat','lon'])

for station in list_of_stations:
    station_locs.loc[station].lat = locs.loc[station].lat
    station_locs.loc[station].lon = locs.loc[station].lon

#Get a new background map figure
fig, ax = plotMap()

sc = ax.scatter(station_locs.lon, station_locs.lat,s=50, marker='o', c='red',
                edgecolors='k',linewidths=0.5, zorder=10, 
                transform=ccrs.PlateCarree())

transform = ccrs.PlateCarree()._as_mpl_transform(ax)
for i in list_of_stations:
    print(list_of_stations[i])
    ax.annotate(list_of_stations[i].split(',')[0],(station_locs.lon[i], station_locs.lat[i]),
                xycoords=transform, transform=ccrs.PlateCarree(),annotation_clip=True, 
                fontsize=12, ha='center',va='center',
                xytext=offsets[i], textcoords='offset points',)

figurePath = '/Users/rantanem/Documents/python/figures/'
figureName = 'arclim_station_map.png'  
plt.savefig(figurePath + figureName,dpi=200,bbox_inches='tight')
   