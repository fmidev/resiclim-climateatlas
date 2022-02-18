#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 14:48:37 2022

This script contains routines which are used to read data from Allas file 
system. 

@author: mprantan
"""
import xarray as xr
import pandas as pd
import allas_utils as allas
import warnings
import sys



def read_daily_data_from_allas(fs, years, yeartype,  var, func):

    
    variables = {'2m_temperature':'t2m',
                 '2m_dewpoint_temperature':'d2m',
                 'snow_cover':'snowc',
                 'total_precipitation':'tp',
                 'snowfall':'sf',
                 'skin_temperature':'skt'}
    
    ## all files
    files = fs.glob('resiclim/daily/'+var+'_'+func+'*')

    da_list = []
    # loop over the years
    for year in years[:]:
        
        # if we want annual (calendar year) values
        if yeartype=='summer':
            matching = [s for s in files if '_'+str(year) in s]
            yearstring = str(year)
        # or from summer to next summer
        elif yeartype=='winter':
            dates = pd.date_range(start=str(year-1)+'-07-01', end=str(year)+'-06-30', freq='MS')
            dates_string = list(dates.strftime('%Y%m'))
            matching = [s for s in files if any(xs in s for xs in dates_string)]
            yearstring = str(year-1)+'-'+str(year)
            
            
        print('Reading ERA5-Land '+func+' '+var+' for '+yearstring, flush=True)
        ds = allas.read_mf_s3(fs, matching, chunks={'time':10}, combine='by_coords')
        # ds = allas.read_mf_s3(fs, matching, chunks={'latitude':90}, combine='by_coords')
    
        da = ds[variables[var]]
        
        da_list.append(da)
    

    da_all = xr.concat(da_list, dim='time')
            
    return da_all

def read_annual_data_from_allas(fs, year, var, func):

    
    variables = {'2m_temperature':'t2m',
                 '2m_dewpoint_temperature':'d2m',
                 'snow_cover':'snowc',
                 'total_precipitation':'tp',
                 'snowfall':'sf'}
    
    ## all files
    files = allas.ls_s3(fs, 'resiclim/coarse/'+var)
    files = fs.glob('resiclim/daily/'+var+'_'+func+'*')

    # select only year 2020
    matching = [s for s in files if str(year) in s]
    
    # open the datasets
    da_list = []
    for f in matching:
       
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            ## if rain, calculate daily sum
            if var=='total_precipitation' or var=='snowfall':
                da = allas.read_mf_s3(fs, f, chunks={'time':250})#[variables[var]].resample(time="1D").sum(dim='time')
            else:
                da = allas.read_mf_s3(fs, f,chunks={'time':250})#[variables[var]].resample(time="1D").mean(dim='time')
            print('Reading '+var+' for',pd.to_datetime(da.time[0].values).strftime('%B %Y'), da.shape)
            da_list.append(da)
        
    da_all = xr.concat(da_list, dim='time')
    
    print('All '+ var + ' fields read!\n')
        
    return da_all

def read_test_data_from_lustre(var, func):
    
    print('Reading ERA5-Land '+var+' from /scratch/')
    
    variables = {'2m_temperature':'t2m',
                 '2m_dewpoint_temperature':'d2m',
                 'snow_cover':'snowc',
                 'total_precipitation':'tp',
                 'snowfall':'sf',
                 'skin_temperature':'skt'}
    
    ## all files
    file = '/scratch/project_2005030/test/era5-land-test_'+variables[var]+'.nc'   

    # open the datasets
    with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            
            if func=='mean':
                da = xr.open_dataset(file)[variables[var]].resample(time="1D").mean(dim='time')
            elif func=='max':
                da = xr.open_dataset(file)[variables[var]].resample(time="1D").max(dim='time')
            elif func=='min':
                da = xr.open_dataset(file)[variables[var]].resample(time="1D").min(dim='time')
            elif func=='sum':
                da = xr.open_dataset(file)[variables[var]].resample(time="1D").sum(dim='time')
            else:
                print('Function not identified')
                sys.exit()
                

        
    return da


    
    

