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
    gs_end = cumsum_da.sel(time=slice(str(y)+'-07-01', str(y)+'-12-31')).argmax(dim='time')
    # Add the missing days from Jan-June
    ndays = len(da_annual.time.sel(time=slice(str(y)+'-01-01', str(y)+'-06-30')))
    gs_end = gs_end + ndays
        

    # calculate the length of GS for each grid point
    gsl = (gs_end - gs_beg).rename('gsl')
        
    # set GSL length zero in areas where the GS did not start before 30th June
    gsl = gsl.where(gs_beg < ndays, 0)

    gsl = gsl.where(gsl>0, 0)*ls_mask.astype(float).assign_coords(time=y) 
    
            
    # Assign attributes
    gsl.attrs['units'] = 'day'
    gsl.attrs['long_name'] = 'Thermal growing season length'
    gsl.attrs['short_name'] = 'GSL'

        
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
    
    # 1/np.nan field (land sea mask)
    ls_mask = da_annual.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
    
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
    gdd = gdd.where(gdd>0, 0)*ls_mask.assign_coords(time=y).rename('gdd').astype(float)
        
    # Assign attributes
    gdd.attrs['units'] = 'C day'
    gdd.attrs['long_name'] = 'Thermal growing degree day sum'
    gdd.attrs['short_name'] = 'GDD'
        
    return gdd

def rain_on_snow(da_tp_winter, da_sf, da_snowc, rain_threshold):
    
    def is_ndjf(month):
        return (month >= 11) | (month <= 2)
    
    def is_mamj(month):
        return (month >= 3) | (month <= 6)
    
    def is_ndjfma(month):
        return (month >= 11) | (month <= 4)
    
    import warnings
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # year
    y = da_tp_winter.time.dt.year[-1].values
    
    print('Calculating the number of rain-on-snow events for '+str(y))        
               
    # Obtain liquid precipitation by subtracting the snowfall from total precipitation
    da_lp_annual = da_tp_winter - da_sf
        
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
    
    # Select only Nov-Apr period of the year
    ros_events = ros_events.sel(time=is_ndjfma(ros_events['time.month']))
        
    # Calculate sum
    ros = ros_events.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    ros = ros.assign_coords(time=y).rename('ros').astype(float)
        
    # Assign attributes
    ros.attrs['units'] = 'events per year'
    ros.attrs['long_name'] = 'Number of rain-on-snow events'
    ros.attrs['short_name'] = 'ROS'
    
    return ros

def rain_on_snow_intensity(da_tp, da_sf, da_snowc, rain_threshold):
    
    def is_ndjfm(month):
        return (month >= 11) | (month <= 3)
    
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
    ros_intensity = ros_intensity.sel(time=is_ndjfm(ros_intensity['time.month']))
    
    # the metric is the cumulative sum of the whole winter year
    rsi = ros_intensity.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    rsi = rsi.assign_coords(time=y).rename('rsi').astype(float)
        
    # Assign attributes
    rsi.attrs['units'] = 'mm days'
    rsi.attrs['long_name'] = 'Total intensity of rain-on-snow events'
    
    return rsi

def winter_warming_events(da_t2mean_winter, da_snowc):
    
    # This function calculates the number of winter warming events. 
    # the events are defined as days in Dec-Mar period, when the grid cell is 
    # snow covered and daily mean temperature rises over 2C. 
    
    def is_ndjfma(month):
        return (month >= 11) | (month <= 4)
    
    import warnings
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
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
    over_two_degrees = da_2t_annual.where(da_2t_annual >= 2.0, np.nan).notnull()
        
    # WW events
    ww_events = (snow_covered * over_two_degrees)
        
    # select only DJFM period
    ww_events = ww_events.sel(time=is_ndjfma(ww_events['time.month']))
    
    # the metric is the cumulative sum over the whole DFJM period
    wwe = ww_events.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    wwe = wwe.assign_coords(time=y).rename('wwint').astype(float)
        
    # Assign attributes
    wwe.attrs['units'] = 'events per year'
    wwe.attrs['long_name'] = 'Number of winter warming events'
    wwe.attrs['short_name'] = 'WWE'
        
    return wwe

def winter_warming_intensity(da_t2mean_winter, da_snowc):
    
    # This function calculates the total intensity of winter warming events. 
    # the events are defined as days in Dec-Mar period, when the grid cell is 
    # snow covered and daily mean temperature rises over 2C. 
    # The intensity of the events is linearly weighted by duration throughout 
    # the event. E.g. for a 3 day event with daily mean air temperatures of 
    # 4 °C, 6 °C and 3 °C, Intensity = (4 ∗ 1) + (6 ∗ 2) + (3 ∗ 3) = 25”.
    # The total intensity is the total cumulative number of all events within
    # a year (Dec-Mar period).
    
    def is_ndjfma(month):
        return (month >= 11) | (month <= 4)
        
    # year
    y = da_t2mean_winter.time.dt.year[-1].values
        
    print('Calculating the intensity of winter warming events for '+str(y)) 
        
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
    ww_intensity = ww_intensity.sel(time=is_ndjfma(ww_intensity['time.month']))
    
    # the metric is the cumulative sum over the whole DFJM period
    ww_accumulative_int = ww_intensity.sum(dim='time', skipna=False)
        
    # Assign coordinate and rename
    wwe_int = ww_accumulative_int.assign_coords(time=y).rename('wwint').astype(float)
        
    # Assign attributes
    wwe_int.attrs['units'] = 'C day'
    wwe_int.attrs['long_name'] = 'Intensity of winter warming events'
    wwe_int.attrs['short_name'] = 'WWI'
        
    return wwe_int

def frost_during_growing_season(da_t2mean_summer, da_skt, basevalue):
    
    # this function calculates the total cumulative frost days during 
    # the growing season. The growing season is determined using so called 
    # integral method (Ruosteenoja et al 2016). The frost is based on 2m-temperature
    
    ## Ruosteenoja et al. (2016): https://doi.org/10.1002/joc.4535
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
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
    
    # FGS is considered positive
    fgs *= -1
        
    # Assign coordinate and rename
    fgs = fgs.assign_coords(time=y).rename('fgs').astype(float)
        
    # Assign attributes
    fgs.attrs['units'] = 'C day'
    fgs.attrs['long_name'] = 'Frost during the growing season'
    fgs.attrs['short_name'] = 'FGS'
        
    return fgs

def vapour_pressure_deficit(da_t2mean_summer, da_d2mean, basevalue):
    
    # This function calculates vapor pressure deficit from 2m temperature 
    # and 2m dew point temperature. Currently the annual 
    # value is obtained by averaging over the whole year.
    #
    
    import warnings
    import math
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
          
    print('Calculating vapour pressure deficit for '+str(y))
        
    # Select annual temperature
    da_2t_annual = da_t2mean_summer - 273.15
    da_2d_annual = da_d2mean - 273.15

    # Calculate Saturated Vapour Pressure in Pa using improved Magnus formula
    # VPsat = (610.7 * 10**((7.5*da_2t_annual)/(237.3+da_2t_annual))) / (1000)
    VPsat = 610.94 * math.exp((17.625*da_2t_annual)/(da_2t_annual + 243.04))
    
    # Calculate actual Vapour Pressure in Pa
    # VPair = (610.7 * 10**((7.5*da_2d_annual)/(237.3+da_2d_annual))) / (1000)
    VPair = 610.94 * math.exp((17.625*da_2d_annual)/(da_2d_annual + 243.04))
        
    # Calculate the deficit
    vpd = VPsat - VPair
    
    #### CALCULATE growing season
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
    # Select VPD within the GS 
    selected_data = vpd * selector # Multidimensional boolean indexing is not supported...
       
    # Calculate mean over the growing season
    vpd = selected_data.where(selected_data > 0).mean(dim='time')
        
    # Assign coordinate and rename
    vpd = vpd.assign_coords(time=y).rename('vpd').astype(float)
    
    # scale to obtain vpd in millibars
    vpd = vpd*10
        
    # Assign attributes
    vpd.attrs['units'] = 'mb'
    vpd.attrs['long_name'] = 'Vapour pressure deficit'
    
    return vpd

def vpd_magnitude_index(da_t2mean_summer, da_d2mean, T90p_vpd, p75vpd, p25vpd):
    
    # This function calculates VPD magnitude index from 2m temperature 
    # and 2m dew point temperature. 
    
    import warnings
    from scipy import ndimage
    
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # take only summer months
    def is_jja(month):
        return (month >= 6) & (month <= 8)
    
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
          
    print('Calculating VPD magnitude index for '+str(y))
    
     # 1/np.nan field (land sea mask)
    ls_mask = da_t2mean_summer.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
        
    # Select annual temperature
    da_2t_annual = da_t2mean_summer - 273.15
    da_2d_annual = da_d2mean - 273.15

    # Calculate Saturated Vapour Pressure in Pa using improved Magnus formula
    VPsat = 610.94 * np.exp((17.625*da_2t_annual)/(da_2t_annual + 243.04))
    
    # Calculate actual Vapour Pressure in Pa
    VPair = 610.94 * np.exp((17.625*da_2d_annual)/(da_2d_annual + 243.04))
        
    # Calculate the deficit
    vpd = VPsat - VPair
    
    # convert threshold coordinates  
    newcoords = pd.to_datetime(y * 1000 + T90p_vpd['doy'], format='%Y%j')   
    T90p_renamed = T90p_vpd.rename({'doy':'time'}).assign_coords(time=newcoords)
    
    # Identify high VPD days
    vpds = (vpd > T90p_renamed) * ls_mask
    
    # generate the structure to label each VPD event
    struct = np.zeros(shape=(3,3,3))
    struct[:, 1, 1] = 1
        
    # label each VPD event
    labels, nb = ndimage.label(vpds, structure=struct)
        
    # calculate the length of each VPD event
    vpd_lengths = np.array(ndimage.sum(vpds, labels, np.arange(labels.max()+1)))
        
    # mask heatwaves which are shorther than three days
    mask = vpd_lengths > 2
    remove_small_vpds = mask[labels.ravel()].reshape(labels.shape)
        
    # make labeled array
    vpd_events = vpd.copy(data=remove_small_vpds)
    
    # select only JJA period
    vpd_values =  vpd.where(vpd_events).sel(time=is_jja(vpd['time.month']))
        
    # calculate the vpd magnitude (based on heatwave magnitude in Dobricic et al. 2020)
    Md = (vpd_values - p25vpd) / (p75vpd - p25vpd)
         
    # calculate sum of the daily magnitudes of Md from the consecutive days composing a vpd event
    # negative Md indices are considered zero
    vpdvalues = (Md.where(Md>0, 0)).values
    cums = np.cumsum(vpdvalues, axis=0)
    weights_by_duration_array = cums - np.maximum.accumulate(cums * (vpdvalues==0), axis=0)  
       
    # make labeled xarray
    cumulative_vpd_magnitude = vpd_values.copy(data=weights_by_duration_array)
        
    # similar to HWM, vpd magnitude index is the maximum value of 
    # Mhw occurring within a given summer 
    vmi = cumulative_vpd_magnitude.max(dim='time', skipna=True)
          
    # Assign coordinate and rename
    vmi = (vmi*ls_mask).assign_coords(time=y).rename('vpd').astype(float)
        
    # Assign attributes
    vmi.attrs['units'] = ''
    vmi.attrs['long_name'] = 'Vapor pressure deficit magnitude index'
    vmi.attrs['short_name'] = 'VPDI'
    
    return vmi

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
    da_t2max_annual = (da_t2max - 273.15)#.sel(latitude=67.6, longitude=133.4).compute()
        
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
    hwi.attrs['units'] = ''
    hwi.attrs['long_name'] = 'Heatwave magnitude index'
    hwi.attrs['short_name'] = 'HWMI'
        
    return hwi
    

def freezing_degree_days(da_t2mean_winter):
    
    ## This function calculates the freezing degree days, using the so-called integral 
    ## method (see Ruosteenoja et al. 2016) to determine the onset and end of the 
    ## freezing season. The freezing season starts when the maximum of cumulative sum of 
    ## daily average temperatures are reached, but not later than Feb 1st. The end
    ## of season is defined when the cumulative minimum of freezing season occurs.
    ## 
    ## The degree days are positive
    ##
    
    ## Ruosteenoja et al. (2016): https://doi.org/10.1002/joc.4535
    
    import warnings 
    # year
    y = da_t2mean_winter.time.dt.year[-1].values
    
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    print('Calculating freezing degree days for '+str(y))
        
    # Select annual temperature and change Kelvin to Celsius
    da_t2mean_annual = da_t2mean_winter - 273.15
    
    # 1/np.nan field (land sea mask)
    ls_mask = da_t2mean_annual.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
    
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
        
    # Define the end of freezing season
    fs_end = cumsum_da.sel(time=slice(str(y)+'-02-01', str(y)+'-06-30')).argmin(dim='time')
    # Select only those locations where the end is non-zero
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
       
    # Calculate sum and multiply by -1 to get positive values
    fdd = selected_data.sum(dim='time', skipna=True) *-1
    
    # Assign coordinate and rename
    fdd = fdd.where(fdd>0, 0)*ls_mask.assign_coords(time=y).rename('fdd').astype(float)
        
    # Assign attributes
    fdd.attrs['units'] = 'C days'
    fdd.attrs['long_name'] = 'Freezing degree days'
    fdd.attrs['short_name'] = 'FDD'
        
    return fdd


def snow_season_length(da_snowc,):
    
    ## This function calculates the length of longest continuous snow season,
    ## defined by the continuous time period between first and last snow date.
    ## The first and last snow date are defined as the first and last days when snow 
    ## cover fraction is > 50 %
    ## The length is given in days
        
    import warnings
    import datetime
    import calendar

    
    def days_in_year(year=datetime.datetime.now().year):
        return 365 + calendar.isleap(year)

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
    ssl = cumulative_snow_length.max(dim='time', skipna=True)
             
    # Assign coordinate and rename
    ssl = (ssl*ls_mask).assign_coords(time=y).rename('ssl').astype(float)

        
    # Assign attributes
    ssl.attrs['units'] = 'day'
    ssl.attrs['long_name'] = 'Snow season length'
    ssl.attrs['short_name'] = 'SSL'

    return ssl

def onset_end_snow_season(da_snowc,):
    
    ## This function calculates the onset and end of continuous snow season,
    ## defined by the continuous time period between first and last snow date.
    ## The first and last snow date are defined as the first and last days when snow 
    ## cover fraction is > 50 %
    ## The length is given in days
        
    import warnings
    import datetime
    import calendar

    
    def days_in_year(year=datetime.datetime.now().year):
        return 365 + calendar.isleap(year)

    # year
    y = da_snowc.time.dt.year[-1].values
    
    # number of days in the previous year
    ndays = days_in_year(y-1)

    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    print('Calculating the onset and end of snow season for '+str(y), flush=True)
           
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
    ssl = cumulative_snow_length.max(dim='time', skipna=True)
    
    # Date of the latest snow-covered day = end of snow season
    ssl_end_date = cumulative_snow_length.idxmax(dim='time', skipna=True)
    ssl_end = ssl_end_date.dt.dayofyear
    
    # Date of the first snow-covered day = onset of snow season
    ssl_beg_date = ssl_end_date - ssl.astype("timedelta64[D]")
    
    # add 365 to the values below 180 
    cond = ssl_beg_date.dt.dayofyear<180
    ssl_beg = xr.where(cond, ssl_beg_date.dt.dayofyear+ndays, ssl_beg_date.dt.dayofyear)
    
    # subtract 365 from the values above 200 
    cond = ssl_end_date.dt.dayofyear>200
    ssl_end = xr.where(cond, ssl_end-ndays, ssl_end)

    # Replace SSL onset and end with NaN on those locations where SSL is 0
    ssl_beg = xr.where(ssl==0, np.nan, ssl_beg)
    ssl_end = xr.where(ssl==0, np.nan, ssl_end)
          
    # Assign coordinate and rename
    ssl_beg = (ssl_beg*ls_mask).assign_coords(time=y).rename('ssb').astype(float)
    ssl_end = (ssl_end*ls_mask).assign_coords(time=y).rename('sse').astype(float)
        
    # Assign attributes
    
    ssl_beg.attrs['units'] = 'day of year'
    ssl_beg.attrs['long_name'] = 'Onset of snow season'
    ssl_beg.attrs['short_name'] = 'SSO'
    
    ssl_end.attrs['units'] = 'day of year'
    ssl_end.attrs['long_name'] = 'End of snow season'
    ssl_end.attrs['short_name'] = 'SSE'

    return ssl_beg, ssl_end


def average_wind_speed(da_u10, da_v10):
    
    ## This function calculates the annual average wind speed 
    
    import warnings 
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    # year
    y = da_u10.time.dt.year[-1].values
    
    print('Calculating the average wind speed for '+str(y), flush=True)

    # calculate wind speed for each day
    ws = np.sqrt(da_u10**2 + da_v10**2)
    
    # calculate annual average
    aws = ws.mean(dim='time').assign_coords(time=y).rename('aws').astype(float)
    
    # Assign attributes
    aws.attrs['units'] = 'm s^-1'
    aws.attrs['long_name'] = 'Annual mean 10-m wind speed'
    aws.attrs['short_name'] = 'WSA'
    
    return aws

def high_wind_events(da_u10max, da_v10max):
    
    ## This function calculates the annual average wind speed 
    
    import warnings 
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
    
    # year
    y = da_u10max.time.dt.year[-1].values
    
    print('Calculating the number of high wind speed events for '+str(y), flush=True)
    
    threshold_ds = xr.open_dataset('/projappl/project_2005030/climateatlas/ws_threshold.nc')
    T90p = threshold_ds.p90
    
    # convert threshold coordinates  
    newcoords = pd.to_datetime(y * 1000 + T90p['doy'], format='%Y%j')   
    T90p_renamed = T90p.rename({'doy':'time'}).assign_coords(time=newcoords)
        
          
    # 1/np.nan field (land sea mask)
    ls_mask = da_u10max.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)

    # calculate wind speed for each day
    ws = np.sqrt(da_u10max**2 + da_v10max**2)
    
    # Identify high WS days
    ws_events = (ws > T90p_renamed) * ls_mask
    
    # Mark grid cells with wind speed < gale with 0
    # ws_events = ws.where((ws > 17.0) | (ws.isnull()), 0)
        
    # Mark grid cells with wind speed > gale with 1
    # ws_events = ws_events.where((ws_events < 17.0) | (ws_events.isnull()),1)
    
    hwe = ws_events.sum(dim='time')*ls_mask.assign_coords(time=y).rename('hwe').astype(float)
    
    # Assign attributes
    hwe.attrs['units'] = 'events per year'
    hwe.attrs['long_name'] = 'Number of high wind speed events'
    hwe.attrs['short_name'] = 'HWE'
    
    return hwe

def annual_mean_temperature(da_t2mean_summer):
    
    ## This function calculates the annual average wind speed 
    
    import warnings 
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
    
    print('Calculating the annual mean temperature for '+str(y), flush=True)

    # calculate annual average
    tavg = da_t2mean_summer.mean(dim='time').assign_coords(time=y).rename('tavg').astype(float)
    
    # Assign attributes
    tavg.attrs['units'] = 'K'
    tavg.attrs['long_name'] = 'Annual mean temperature'
    tavg.attrs['short_name'] = 'TAVG'
    tavg.attrs['standard_name'] = 'air_temperature'
    
    return tavg

def annual_precipitation(da_tp_summer):
    
    ## This function calculates the annual precipitation 
    
    import warnings 
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    # year
    y = da_tp_summer.time.dt.year[-1].values
    
    print('Calculating the annual precipitation for '+str(y), flush=True)

    # calculate annual average
    tpa = da_tp_summer.sum(dim='time', skipna=False).assign_coords(time=y).rename('tpa').astype(float)
    
    # Assign attributes
    tpa.attrs['units'] = 'm'
    tpa.attrs['long_name'] = 'Annual precipitation'
    tpa.attrs['short_name'] = 'PRA'
    
    return tpa

def annual_snowfall(da_sf):
    
    ## This function calculates the annual average wind speed 
    
    import warnings 
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    # year
    y = da_sf.time.dt.year[-1].values
    
    print('Calculating the annual snowfall for '+str(y), flush=True)

    # calculate annual average
    sfa = da_sf.sum(dim='time', skipna=False).assign_coords(time=y).rename('sfa').astype(float)
    
    # Assign attributes
    sfa.attrs['units'] = 'm'
    sfa.attrs['long_name'] = 'Annual snowfall'
    sfa.attrs['short_name'] = 'SFA'
    
    return sfa

def summer_warmth_index(da_t2mean_summer):
    
    ## This function calculates the annual average wind speed 
    
    import warnings 
    
    warnings.simplefilter("ignore", category=RuntimeWarning)
        
    # year
    y = da_t2mean_summer.time.dt.year[-1].values
    
    print('Calculating the summer warmth index for '+str(y), flush=True)

    # calculate monthly averages
    t2mean_monthly = da_t2mean_summer.groupby(da_t2mean_summer.time.dt.month).mean() -273.15
    
    # months which are above 0C
    summer_months = t2mean_monthly.where(t2mean_monthly >= 0.0)
    
    # 1/np.nan field (land sea mask)
    ls_mask = da_t2mean_summer.isel(time=0).notnull()
    ls_mask = ls_mask.where(ls_mask, np.nan)
    
    # calculate the index
    swi = summer_months.sum(dim='month')*ls_mask
    
    swi = swi.assign_coords(time=y).rename('swi').astype(float)
    
    # Assign attributes
    swi.attrs['units'] = 'C'
    swi.attrs['long_name'] = 'Summer warmth index'
    swi.attrs['short_name'] = 'SWI'
    
    return swi