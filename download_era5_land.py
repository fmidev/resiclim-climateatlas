#!/usr/bin/env python
#
#

import os, sys, itertools
import numpy as np
import cdsapi
server = cdsapi.Client()

import allas_utils as allas


varname = str(sys.argv[1]) 


years = np.arange(1950,2022).astype(int)
years = np.arange(2021,1949,-1).astype(int)
#years = (2020,)
months = np.arange(1,13)

for month, year in itertools.product(months,years):
    
    basename = '%s_era5Land_%04d%02d' % (varname,year,month)
    nc_file  = '%s.nc'  % (basename)
    remote_path = 'resiclim/'+varname+'/'+nc_file
    
    if allas.isfile_s3(remote_path):
        pass
    
    else:
        
        opts = {
                'product_type'  : 'reanalysis',  
                'area'          : [90, -180, 45, 180,],
                'variable'      : varname, 
                'year'          : '%04d' % (year),
                'month'         : '%02d' % (month),
                'day'           : ['%02d' % (i+1) for i in range(31)], 
                'time'          : ['%02d:00' % (i) for i in range(24)], 
                'format'        : 'netcdf',
               }
        
        print('Fetching data for', nc_file)
        
        # Retrieve from Climate Data Store
        try:
            server.retrieve('reanalysis-era5-land', opts, nc_file)
        except:
            print('Retrieval failed for',nc_file)
            
        
        # Upload to S3 Allas and remove the temporal file from disk
        try:
            allas.upload_data_s3(remote_path, nc_file)
            os.remove(nc_file)
        except:
            print('Upload failed for',nc_file)

