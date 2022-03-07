#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 12:28:43 2022

@author: mprantan
"""
import xarray as xr
import pandas as pd
import  numpy as np
import warnings

def thermal_growing_season_length(da_t2mean_summer, basevalue):
    
    ## This function calculates the length of growing season, using the so-called integral 
    ## method (see Ruosteenoja et al. 2016), which identifies the date after the absolute 
    ## minimum of the sum(Tday-basevalue) has been reached (gs_beg) and analogously 
    ## gs_end when the absolute maximum of the sum(Tday-basevalue) has been reached, 
    ## but not earlier than 1st of July.
    ## Ruosteenoja et al. (2016): https://doi.org/10.1002/joc.4535
    ## The length is given in days
    ## Preferable basevalue in the Arctic is 5 C or 3 C
    
    
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
        
    print('Calculating thermal growing season lenght for '+str(y), flush=True)
        
    # Select annual temperature and convert from Kelvin to Celsius
    da_annual = da_t2mean_summer - 273.15#.where(da_t2mean_summer.time.dt.year == y, drop=True)#.sel(latitude=71.5, longitude=-180)

    # subtract the base value
    da_annual -= basevalue
    
    # 1/np.nan field (land sea mask)
    ls_mask = da_annual.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
    
    # cumulative temperature sum
    cumsum_da = da_annual.cumsum(dim='time', skipna=False)
        
    # fill the nans of the array with obviously wrong value
    cumsum_da = cumsum_da.fillna(-99999)
        
    # day of minimum before June 30th + 1 day = beginning of GS
    day_min = cumsum_da.sel(time=slice(str(y)+'-01-01', str(y)+'-06-30')).argmin(dim='time', skipna=False)
    gs_beg = day_min + 1
        
    # gs_beg = gs_beg.where(gs_beg>1)
        
    # day of maximum = end of GS
    gs_end = cumsum_da.sel(time=slice(str(y)+'-06-01', str(y)+'-12-31')).argmax(dim='time')
    # Add the missing days from Jan-May
    gs_end = gs_end + len(da_annual.time.sel(time=slice(str(y)+'-01-01', str(y)+'-05-31')))
        

    # calculate the length of GS for each grid point
    gsl = (gs_end - gs_beg)*ls_mask.assign_coords(time=y).rename('gsl')
        
    # replace the wrong values with NaN
    gsl = gsl.where(gsl > 0).astype(float) 
        
    # Assign attributes
    gsl.attrs['long_name'] = 'Length of thermal growing season in days'

        
    return gsl

def thermal_growing_degree_days(da_t2mean_summer, basevalue):
    
    ## This function calculates the growing degree days, using the so-called integral 
    ## method (see Ruosteenoja et al. 2016), which identifies the date after the absolute 
    ## minimum of the sum(Tday-basevalue) has been reached (gs_beg) and analogously 
    ## gs_end when the absolute maximum of the sum(Tday-basevalue) has been reached, 
    ## but not earlier than 1st of July.
    ## Ruosteenoja et al. (2016): https://doi.org/10.1002/joc.4535
    
    import warnings
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
        
    print('Calculating thermal growing degree days for '+str(y))
        
    # Select annual temperature in celsius
    da_annual = da_t2mean_summer - 273.15 #.where(da_2t.time.dt.year == y, drop=True)#.sel(latitude=71.5, longitude=-180)
        
    # subtract the base value
    da_annual -= basevalue
    
    # cumulative temperature sum
    cumsum_da = da_annual.cumsum(dim='time', skipna=False)
        
    # fill the nans of the array with obviously wrong value
    cumsum_da = cumsum_da.fillna(-99999)
        
    # lenght of the period from which the minimum is searched (in days)
    # January to the end of June
    len_period = len(da_annual.time.sel(time=slice(str(y)+'-01-01', str(y)+'-06-30')))
        
    # Define the beginning of growing season
    day_min = cumsum_da.sel(time=slice(str(y)+'-01-01', str(y)+'-06-30')).argmin(dim='time')
    gs_beg = day_min + 1
        
    # Define the end of growing season
    gs_end = cumsum_da.sel(time=slice(str(y)+'-07-01', str(y)+'-12-31')).argmax(dim='time')
    # Select only those location where the end is non-zero
    gs_end = gs_end.where(gs_end>0)
        
    # Add the missing days from Jan-May to get the actual day of year
    gs_end = gs_end + len_period
        
    # Create a helper time array - each element's value is the timestamp's value
    time = da_annual.coords['time'].dt.dayofyear
    expanded_time = time.expand_dims({'latitude': da_annual.latitude, 'longitude':da_annual.longitude})
        
    # where() -- for each element, if condition is false, set element to nan
    e1 = expanded_time.where(expanded_time <= gs_end, np.nan)
    e2 = e1.where(e1 >= gs_beg, np.nan)
    # make 1/np.nan array
    selector = e2.where(e2.isnull(),1)  
        
    # Now that we have an indexer array that selects the elements we want, we can calculate our result
    # Negative days within the season do not reduce the sum; thus replace below-zero temperatures by zero
    selected_data = da_annual.where(da_annual>=0,0) * selector # Multidimensional boolean indexing is not supported...
       
    # Calculate sum
    gdd = selected_data.sum(dim='time', skipna=True)
        
    # Assign coordinate and rename
    gdd = gdd.where(gdd>0).assign_coords(time=y).rename('gdd').astype(float)
        
    # Assign attributes
    gdd.attrs['units'] = 'C day'
    gdd.attrs['long_name'] = 'Growing degree day sum'
        
    return gdd

def rain_on_snow(da_tp, da_sf, da_snowc, rain_threshold):
    
    def is_djfm(month):
        return (month >= 12) | (month <= 3)
    
    # this function calculates rain-on-snow events
    
    import warnings
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # year
    y = da_tp.time.dt.year[-1].values
    
    print('Calculating the number of rain-on-snow events for '+str(y))        
               
    # Obtain liquid precipitation by subtracting the snowfall from total precipitation
    da_lp_annual = da_tp - da_sf
        
    # When the snow cover is 0.5 or greater in the grid cell,
    # we consider it snow-covered. Retain nan-points over the sea
        
    # Mark grid cells with snow cover < 50 with 0
    snow_covered = da_snowc.where((da_snowc > 50) | (da_snowc.isnull()), 0)
        
    # Mark grid cells with snow cover > 50 with 1
    snow_covered = snow_covered.where((snow_covered < 50) | (snow_covered.isnull()),1)
        
    # Liquid precipitation needs to be higher than the threshold
    rain = da_lp_annual.where(da_lp_annual > rain_threshold, np.nan).notnull()
        
    # ROS events
    ros_events = (snow_covered * rain)
    
    # Select only DJFM period
    ros_events = ros_events.sel(time=is_djfm(ros_events['time.month']))
        
    # Calculate sum
    ros = ros_events.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    ros = ros.assign_coords(time=y).rename('ros').astype(float)
        
    # Assign attributes
    ros.attrs['units'] = 'events per year'
    ros.attrs['long_name'] = 'Rain-on-snow events'
    
    return ros

def rain_on_snow_intensity(da_tp, da_sf, da_snowc, rain_threshold):
    
    def is_djfm(month):
        return (month >= 12) | (month <= 3)
    
    # this function calculates rain-on-snow events
    
    import warnings
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # year
    y = da_tp.time.dt.year[-1].values
    
    print('Calculating the intensity of rain-on-snow events for '+str(y))        
               
    # Obtain liquid precipitation by subtracting the snowfall from total precipitation
    da_lp_annual = da_tp - da_sf
        
    # When the snow cover is 0.5 or greater in the grid cell,
    # we consider it snow-covered. Retain nan-points over the sea
        
    # Mark grid cells with snow cover < 50 with 0
    snow_covered = da_snowc.where((da_snowc > 50) | (da_snowc.isnull()), 0)
        
    # Mark grid cells with snow cover > 50 with 1
    snow_covered = snow_covered.where((snow_covered < 50) | (snow_covered.isnull()),1)
        
    # Liquid precipitation needs to be higher than the threshold
    rain = da_lp_annual.where(da_lp_annual > rain_threshold, np.nan).notnull()
        
    # ROS events
    ros_events = (snow_covered * rain)
    
    # load all values (NOTE: this may take long!)
    ros_values = ros_events.values
        
    # calculate weights by event duration
    cums = np.cumsum(ros_values, axis=0)
    weights_by_duration_array = cums - np.maximum.accumulate(cums * (ros_values == 0), axis=0)  
    weights_by_duration = ros_events.copy(data=weights_by_duration_array)
        
        
    # Calculate intensity of events
    ros_intensity = (da_lp_annual - rain_threshold) * ros_events * weights_by_duration
        
    # select only DJFM period
    ros_intensity = ros_intensity.sel(time=is_djfm(ros_intensity['time.month']))
    
    # the metric is the cumulative sum of the whole winter year
    rsi = ros_intensity.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    rsi = rsi.assign_coords(time=y).rename('rsi').astype(float)
        
    # Assign attributes
    rsi.attrs['units'] = 'mm days'
    rsi.attrs['long_name'] = 'Total intensity of rain-on-snow events'
    
    return rsi


def winter_warming(da_t2mean_winter, da_snowc):
    
    # This function calculates the total intensity of winter warming events. 
    # the events are defined as days in Dec-Mar period, when the grid cell is 
    # snow covered and daily mean temperature rises over 2C. 
    # The intensity of the events is linearly weighted by duration throughout 
    # the event. E.g. for a 3 day event with daily mean air temperatures of 
    # 4 °C, 6 °C and 3 °C, Intensity = (4 ∗ 1) + (6 ∗ 2) + (3 ∗ 3) = 25”.
    # The total intensity is the total cumulative number of all events within
    # a year (Dec-Mar period).
    
    def is_djfm(month):
        return (month >= 12) | (month <= 3)
        
    # year
    y = da_t2mean_winter.time.dt.year[-1].values
        
    print('Calculating winter warming events for '+str(y)) 
        
    # select annual data and convert from Kelvin to Celsius
    da_2t_annual = da_t2mean_winter - 273.15 

        
    # When the snow cover is 0.5 or greater in the grid cell,
    # we consider it snow-covered. 
        
    # Replace grid cells with snow cover < 50 by 0
    snow_covered = da_snowc.where((da_snowc > 50) | (da_snowc.isnull()), 0)
        
    # Replace grid cells with snow cover > 50 with 1
    snow_covered = snow_covered.where((snow_covered < 50) | (snow_covered.isnull()),1)
        
    # Daily mean temperature needs to be at least 2 C
    over_two_degrees = da_2t_annual.where(da_2t_annual > 2, np.nan).notnull()
        
    # WW events
    ww_events = (snow_covered * over_two_degrees)
        
    # load all values (NOTE: this may take long!)
    warnings.simplefilter("ignore", category=RuntimeWarning)
    ww_values = ww_events.values
        
    # calculate weights by event duration
    cums = np.cumsum(ww_values, axis=0)
    weights_by_duration_array = cums - np.maximum.accumulate(cums * (ww_values == 0), axis=0)  
    weights_by_duration = ww_events.copy(data=weights_by_duration_array)
        
        
    # Calculate intensity of events
    ww_intensity = (da_2t_annual - 2) * ww_events * weights_by_duration
        
    # select only DJFM period
    ww_intensity = ww_intensity.sel(time=is_djfm(ww_intensity['time.month']))
    
    # the metric is the cumulative sum over the whole DFJM period
    ww_accumulative_int = ww_intensity.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    wwe_int = ww_accumulative_int.assign_coords(time=y).rename('wwint').astype(float)
        
    # Assign attributes
    wwe_int.attrs['units'] = 'degrees'
    wwe_int.attrs['long_name'] = 'Total intensity of winter warming events'
        
    return wwe_int

def frost_during_growing_season(da_t2mean_summer, da_skt, basevalue):
    
    # this function calculates the total cumulative frost days during 
    # the growing season. The growing season is determined using so called 
    # integral method (Ruosteenoja et al 2016). The frost is based on 2m-temperature
    
    ## Ruosteenoja et al. (2016): https://doi.org/10.1002/joc.4535
    
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
        
    print('Calculating frost during growing season for '+str(y))
        
    # Select annual temperature
    da_2t_annual = da_t2mean_summer - 273.15 
    da_skt_annual = da_skt - 273.15 
    
    # subtract the base value
    da_2t_gs = da_2t_annual- basevalue
    
    # cumulative temperature sum
    cumsum_da = da_2t_gs.cumsum(dim='time', skipna=False)
        
    # fill the nans of the array with obviously wrong value
    cumsum_da = cumsum_da.fillna(-99999)
        
    # lenght of the period from which the minimum is searched (in days)
    # January to the end of June
    len_period = len(da_2t_annual.time.sel(time=slice(str(y)+'-01-01', str(y)+'-06-30')))
        
    # Define the beginning of growing season
    day_min = cumsum_da.sel(time=slice(str(y)+'-01-01', str(y)+'-06-30')).argmin(dim='time')
    gs_beg = day_min + 1
        
    # Define the end of growing season
    gs_end = cumsum_da.sel(time=slice(str(y)+'-07-01', str(y)+'-12-31')).argmax(dim='time')
    # Select only those location where the end is non-zero
    gs_end = gs_end.where(gs_end>0)
        
    # Add the missing days from Jan-May to get the actual day of year
    gs_end = gs_end + len_period
        
    # Create a helper time array - each element's value is the timestamp's value
    time = da_2t_annual.coords['time'].dt.dayofyear
    expanded_time = time.expand_dims({'latitude': da_2t_annual.latitude, 
                                      'longitude':da_2t_annual.longitude})
        
    # where() -- for each element, if condition is false, set element to nan
    e1 = expanded_time.where(expanded_time <= gs_end, 0.0)
    e2 = e1.where(e1 >= gs_beg, 0.0)
        
    # make 1/np.nan array
    selector = e2.where(e2==0.0,1.)  
        
    # Now that we have an indexer array that selects the elements we want, we can calculate our result
    # Select those days within the GS when temperature is below zero
    selected_data = da_skt_annual.where((da_skt_annual < 0.0) | (da_skt_annual.isnull()), 0.0) * selector # Multidimensional boolean indexing is not supported...
       
    # Calculate sum
    fgs = selected_data.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    fgs = fgs.assign_coords(time=y).rename('fgs').astype(float)
        
    # Assign attributes
    fgs.attrs['units'] = 'C day'
    fgs.attrs['long_name'] = 'Frost during the growing season'
        
    return fgs

def vapour_pressure_deficit(da_t2mean_summer, da_d2mean):
    
    # This function calculates vapor pressure deficit from 2m temperature 
    # and 2m dew point temperature. Currently the annual 
    # value is obtained by averaging over the whole year.
    #
    
    import warnings
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
          
    print('Calculating vapour pressure deficit for '+str(y))
        
    # Select annual temperature
    da_2t_annual = da_t2mean_summer - 273.15
    da_2d_annual = da_d2mean - 273.15

    # Calculate Saturated Vapour Pressure in kPa
    VPsat = (610.7 * 10**((7.5*da_2t_annual)/(237.3+da_2t_annual))) / (1000)
        
    # Calculate actual Vapour Pressure in kPa
    VPair = (610.7 * 10**((7.5*da_2d_annual)/(237.3+da_2d_annual))) / (1000)
        
    # Calculate the deficit
    vpd = VPsat - VPair
        
    # Calculate annual mean
    vpd = vpd.mean(dim='time', skipna=True)
        
    # Assign coordinate and rename
    vpd = vpd.assign_coords(time=y).rename('vpd').astype(float)
        
    # Assign attributes
    vpd.attrs['units'] = 'kPa'
    vpd.attrs['long_name'] = 'Vapour pressure deficit'
    
    return vpd

def heatwave_magnitude_index(da_t2max, T90p, p75max, p25max):
    
    # This function calculates the annual heatwave magnitude index which is 
    # described in Dobricic et al (2020) and Russo et al (2015).
    
    # Dobricic et al (2020): https://iopscience.iop.org/article/10.1088/1748-9326/ab6398/meta
    # Russo et al (2015): https://iopscience.iop.org/article/10.1088/1748-9326/10/12/124003
    
    from scipy import ndimage
    import warnings
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # take only summer months
    def is_jja(month):
        return (month >= 6) & (month <= 8)
    
    # year
    y = da_t2max.time.dt.year[-1].values
        
    print('Calculating heatwave magnitude index for '+str(y))
        
    # Select annual temperature
    da_t2max_annual = da_t2max - 273.15 #.where(da_t2max.time.dt.year == y, drop=True)#.sel(latitude=61, longitude=27)
        
    # 1/np.nan field (land sea mask)
    ls_mask = da_t2max_annual.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
        
    # convert threshold coordinates  
    newcoords = pd.to_datetime(y * 1000 + T90p['doy'], format='%Y%j')   
    T90p_renamed = T90p.rename({'doy':'time'}).assign_coords(time=newcoords)
        
    # Identify heatwave days
    heatwaves = (da_t2max_annual > T90p_renamed) * ls_mask
        
    # generate the structure to label each heatwave event
    struct = np.zeros(shape=(3,3,3))
    struct[:, 1, 1] = 1
        
    # label each heatwave event
    labels, nb = ndimage.label(heatwaves, structure=struct)
        
    # calculate the length of each heatwave
    heatwave_lengths = np.array(ndimage.sum(heatwaves, labels, np.arange(labels.max()+1)))
        
    # mask heatwaves which are shorther than three days
    mask = heatwave_lengths > 2
    remove_small_heatwaves = mask[labels.ravel()].reshape(labels.shape)
        
    # make labeled array
    heatwave_events = da_t2max_annual.copy(data=remove_small_heatwaves)
        
    # select only JJA period
    heatwave_temps =  da_t2max_annual.where(heatwave_events).sel(time=is_jja(da_t2max_annual['time.month']))
        
    # calculate the heatwave magnitude (based on Dobricic et al. 2020)
    Md = (heatwave_temps - p25max) / (p75max - p25max)
         
    # calculate sum of the daily magnitudes of Md from the consecutive days composing a heat wave
    # negative Md indices are considered zero
    heatwavevalues = (Md.where(Md>0, 0)).values
    cums = np.cumsum(heatwavevalues, axis=0)
    weights_by_duration_array = cums - np.maximum.accumulate(cums * (heatwavevalues==0), axis=0)  
       
    # make labeled xarray
    cumulative_heatwave_magnitude = heatwave_temps.copy(data=weights_by_duration_array)
        
    # Based on Dobricic et al. 2020, heatwave index is the maximum value of 
    # Mhw occurring within a given summer 
    hwi = cumulative_heatwave_magnitude.max(dim='time', skipna=True)
          
    # Assign coordinate and rename
    hwi = (hwi*ls_mask).assign_coords(time=y).rename('hwi').astype(float)
        
    # Assign attributes
    hwi.attrs['units'] = ' '
    hwi.attrs['long_name'] = 'Heatwave magnitude index'
        
    return hwi
    

def freezing_degree_days(da_t2mean_winter):
    
    ## This function calculates the freezing degree days, using the so-called integral 
    ## method (see Ruosteenoja et al. 2016) to determine the onset and end of the 
    ## freezing season. The freezing season starts when the maximum of cumulative sum of 
    ## daily average temperatures are reached, but not later than Feb 1st. The end
    ## of season is defined when the cumulative minimum of freezing season occurs.
    ## 
    ## The degree days are negative
    ##
    
    ## Ruosteenoja et al. (2016): https://doi.org/10.1002/joc.4535
    
    import warnings 
    # year
    y = da_t2mean_winter.time.dt.year[-1].values
    
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    print('Calculating freezing degree days for '+str(y))
        
    # Select annual temperature and change Kelvin to Celsius
    da_t2mean_annual = da_t2mean_winter - 273.15
    
     # cumulative temperature sum
    cumsum_da = da_t2mean_annual.cumsum(dim='time', skipna=False)
        
    # fill the nans of the array with obviously wrong value
    cumsum_da = cumsum_da.fillna(-99999)
        
    # lenght of the period from which the maximum is searched (in days)
    # July to the end of December
    len_period = len(da_t2mean_annual.time.sel(time=slice(str(y-1)+'-07-01', str(y)+'-01-31')))
        
    # Define the beginning of freezing season
    day_max = cumsum_da.sel(time=slice(str(y-1)+'-07-01', str(y)+'-01-31')).argmax(dim='time')
    fs_beg = day_max + 1
        
    # Define the end of growing season
    fs_end = cumsum_da.sel(time=slice(str(y)+'-02-01', str(y)+'-06-30')).argmin(dim='time')
    # Select only those location where the end is non-zero
    fs_end = fs_end.where(fs_end>0)
        
    # Add the missing days from Jul-Jan to get the actual day of year
    fs_end = fs_end + len_period
    
    # Create a helper time array - each element's value is the timestamp's value    
    time = da_t2mean_annual.coords['time'].copy(data=np.arange(1, np.shape(da_t2mean_annual)[0]+1))
    expanded_time = time.expand_dims({'latitude': da_t2mean_annual.latitude, 
                                      'longitude':da_t2mean_annual.longitude})
        
    # where() -- for each element, if condition is false, set element to nan
    e1 = expanded_time.where(expanded_time <= fs_end, np.nan)
    e2 = e1.where(e1 >= fs_beg, np.nan)
    # make 1/np.nan array
    selector = e2.where(e2.isnull(),1)  
        
    # Now that we have an indexer array that selects the elements we want, we can calculate our result
    # Positive days within the freezing season do not reduce the sum; 
    # thus replace above-zero temperatures by zero
    selected_data = da_t2mean_annual.where(da_t2mean_annual<=0,0) * selector # Multidimensional boolean indexing is not supported...
       
    # Calculate sum
    fdd = selected_data.sum(dim='time', skipna=True)
    
    # Assign coordinate and rename
    fdd = fdd.where(fdd<0).assign_coords(time=y).rename('fdd').astype(float)
        
    # Assign attributes
    fdd.attrs['units'] = 'C day'
    fdd.attrs['long_name'] = 'Freezing degree day sum'
        
    return fdd

def snow_season_length(da_snowc, ):
    
    ## This function calculates the length of snow season,
    ## defined by the time period between first and last snow date.
    ## The first and last snow date are defined as the first and last days when snow 
    ## cover fraction is > 50 %
    ## The length is given in days
    
    # year
    y = da_snowc.time.dt.year[-1].values
        
    print('Calculating snow season lenght for '+str(y), flush=True)
    
    # 1/np.nan field (land sea mask)
    ls_mask = da_snowc.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
        
    # Mark grid cells with snow cover < 50 with 0
    snow_covered = da_snowc.where((da_snowc > 50) | (da_snowc.isnull()), 0)
        
    # Mark grid cells with snow cover > 50 with 1
    snow_covered = snow_covered.where((snow_covered < 50) | (snow_covered.isnull()),1)
    
    # fill the nans of the array with obviously wrong value
    snow_covered = snow_covered.fillna(-99999)
    
    # first day of snow
    fds = snow_covered.argmax(dim='time')
    
    # last day of snow
    lds = len(snow_covered.time) - snow_covered.reindex(time=snow_covered.time[::-1]).argmax(dim='time')
    
    # length of snow season
    lss = (lds - fds).assign_coords(time=y).rename('lss')
    
    # replace the wrong values with NaN
    lss = lss.where(lss < 365).astype(float)
        
    # Assign attributes
    lss.attrs['units'] = ' '
    lss.attrs['long_name'] = 'Length of snow season in days'

    return lss

def longest_snow_period(da_snowc,):
    
    ## This function calculates the length of longest continuous snow season,
    ## defined by the continuous time period between first and last snow date.
    ## The first and last snow date are defined as the first and last days when snow 
    ## cover fraction is > 50 %
    ## The length is given in days
        
    import warnings 
    
    # year
    y = da_snowc.time.dt.year[-1].values
    

    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    print('Calculating the longest continuous snow period for '+str(y), flush=True)
           
    # 1/np.nan field (land sea mask)
    ls_mask = da_snowc.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)

    # Mark grid cells with snow cover < 50 with 0
    snow_covered = da_snowc.where((da_snowc > 50) | (da_snowc.isnull()), 0)
        
    # Mark grid cells with snow cover > 50 with 1
    snow_covered = snow_covered.where((snow_covered < 50) | (snow_covered.isnull()),1)
        
    # calculate cumulative sum of snow covered days
    snowvalues = snow_covered.values
    cums = np.cumsum(snowvalues, axis=0)
    weights_by_duration_array = cums - np.maximum.accumulate(cums * (snowvalues==0), axis=0)  
    
    # make labeled xarray
    cumulative_snow_length = snow_covered.copy(data=weights_by_duration_array)
        
    # The longest continuous snow period is the maximum of the cumulative sum
    lsp = cumulative_snow_length.max(dim='time', skipna=True)
          
    # Assign coordinate and rename
    lsp = (lsp*ls_mask).assign_coords(time=y).rename('lcs').astype(float)
        
    # Assign attributes
    lsp.attrs['units'] = ' '
    lsp.attrs['long_name'] = 'The longest continuous period of snow'

    return lsp
