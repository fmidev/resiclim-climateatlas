#!/usr/bin/env python





import matplotlib.pyplot as plt
import allas_utils as allas
import xarray as xr

import os
import s3fs, boto3






# Set up the S3 file system
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 

fs = s3fs.S3FileSystem(anon=False, key=access, secret=secret,
    client_kwargs={'endpoint_url': 'https://a3s.fi'})

'''
s3 = boto3.resource(service_name='s3', 
                    aws_access_key_id=access, 
                    aws_secret_access_key=secret, 
                    endpoint_url='https://a3s.fi')

'''






# List all precipitation files
files = allas.ls_s3(fs, 'resiclim/total_precipitation/')



# select only year 2020
matching = [s for s in files if "2020" in s]

# Read one year of data lazily.
# Chunk size was selected based on this: (365 days * 24 hours) / (40 computing cores) = 220
ds = allas.read_mf_s3(fs, matching, chunks={'time':250})



# Perform a set of operations lazily
ds['tp'] += 5
mean = ds['tp'].mean(dim='time')


# Finally compute everything in parallel
result = mean.compute()




