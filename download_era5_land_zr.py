#!/usr/bin/env python
#
#

import os, sys, itertools
import numpy as np
import xarray as xr

import cdsapi
server = cdsapi.Client()


sys.path.append('/users/kamarain/resiclim-climateatlas')
import allas_utils as allas


varname = str(sys.argv[1]) 


years = np.arange(1950,2022).astype(int)
years = np.arange(2021,1949,-1).astype(int)
years = np.arange(2021,2015,-1).astype(int)
years = (2020,)
months = np.arange(1,13)

for year, month in itertools.product(years, months):
    
    basename = '%s_era5Land_%04d%02d' % (varname,year,month)
    nc_file  = '%s.nc'  % (basename)
    zr_file  = '%s.zarr'  % (basename)
    
    remote_path = 's3://resiclim/test/'+varname+'/'+zr_file
    
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
                #'grid'          : '0.5/0.5', #'0.125/0.125',
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
            ds = xr.open_dataset(nc_file)
            task = allas.xr_to_zarr_s3(ds, remote_path)
            os.remove(nc_file)
        except:
            print('Upload failed for',nc_file)

