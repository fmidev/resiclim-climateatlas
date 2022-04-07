#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:00:19 2022

This script contains routines which read data from FMI Open data.

### WORKS ONLY IN FMI INTERNAL NETWORK!!!  ###

@author: rantanem
"""

#
import requests
import os
from datetime import datetime
import pandas as pd
import sys



def update_station_data(station):
    
    
    # station id's for selected stations
    stationIDs = {'kaisaniemi':100971,
                  'sodankyla':101932, 
                  'ahtari': 101520}
    
    # retrieve daily average, maximum and minimum temperature + precipitation    
    var_longnames = {'TA_P1D_AVG':'Average temperature',
                     'TA_PT24H_MAX':'Maximum temperature',
                     'TA_PT24H_MIN':'Minimum temperature',
                     'PRA_PT24H_SUM':'Precipitation',
                     'SND_P1D_INSTANT':'Snow depth',}

   
    ID = stationIDs[station]
    dataset = None
    start_time = "19600101T0000"
    end_time = "20211231T2359"
    
    stationfilename=str(station) + '_updated.csv'
    df_path = '/Users/rantanem/Documents/python/resiclim-climateatlas/validation/'+stationfilename
    
    # check if the file exists and is not empty
    if os.path.exists(df_path) and os.path.getsize(df_path) > 0:
        # if yes, read the data
        df = pd.read_csv(df_path)
    
    else:
        df, stationname, lat, lon = get_station_data(ID, var_longnames, start_time, end_time)
    
    
    # save csv-file
    df_csv = df.copy()
    df_csv.to_csv('/Users/rantanem/Documents/python/resiclim-climateatlas/validation/data/'+stationfilename, na_rep='NaN',index=False)
    
    df.index = pd.to_datetime(df[['Year','Month','Day']])


    return df


def get_station_data(stationID, variables, start_time, end_time ):
    

    # reading parameters
    # ==============
    asemaFMISID = [stationID]  # lista haettavista asemista (voi olla useita)
    suurekoodit = ["utctime", 'stationname', 'latitude', 'longitude'] + list(variables.keys())  # aikaleima = "utctime"
    alkuaika = start_time  # muodossa YYYYMMDDTHHMM
    loppuaika = end_time  # muodossa YYYYMMDDTHHMM
    tallennuskansio = "data"
    # # #
    
    os.makedirs(tallennuskansio, exist_ok=True)
    

    for asema in asemaFMISID:
        print('  Haetaan {}...'.format(asema), end='', flush=True)
        params = {
            "fmisid": asema,
            "producer": "observations_fmi",
            "param": ",".join(suurekoodit),
            "starttime": alkuaika,
            "endtime": loppuaika,
            "format": "ascii",
            "separator": ",",
            "precision": "double"
            }
        r = requests.get('http://smartmet.fmi.fi/timeseries', params=params)
        if r.ok:
            if len(r.text) == 0:
                print('Ei dataa.')
            else:
                with open(os.path.join(tallennuskansio, "{}.csv".format(asema)), 'wt') as f:
                    f.write("# Data FMI SMARTMET-tietokannasta,Haettu {}Z\n".format(datetime.utcnow().strftime("%d.%m.%Y %H:%M:%S")))
                    f.write("# Asema FMISID: {}\n".format(asema))
                    f.write(",".join(suurekoodit) + '\n')
                    f.write(r.text)
                    f.write('\n')
                    print('Valmis!')
        else:
            sys.exit('Haku epäonnistui.')
            

    ### put data into dataframe
    df = pd.read_csv(os.path.join(tallennuskansio, "{}.csv".format(asema)), header=2)
    stationname = df['stationname'][0]
    lat = df['latitude'][0]
    lon = df['longitude'][0]
    
    df['Year'] = pd.to_datetime(df.utctime.values).year
    df['Month'] = pd.to_datetime(df.utctime.values).month
    df['Day'] = pd.to_datetime(df.utctime.values).day
    
    
    df.drop(columns=['utctime'], inplace=True)
    df = df[['Year', 'Month', 'Day']+list(variables.keys())]
    
    ### mark -1's to zeros in precipitation
    if 'PRA_PT24H_SUM' in df:
        df['PRA_PT24H_SUM'][df['PRA_PT24H_SUM']==-1] = 0
    
    # rename columns
    df.rename(columns=variables, inplace=True)
    
    return df, stationname, lat, lon
