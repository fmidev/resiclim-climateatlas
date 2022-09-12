#!/usr/bin/env python



import io, os
import s3fs, boto3, fsspec
import xarray as xr
import pandas as pd

'''
# Set up the S3 file system 
access = os.getenv('S3_RESICLIM_ACCESS') 
secret = os.getenv('S3_RESICLIM_SECRET') 


# Global handles
FS = s3fs.S3FileSystem(anon=False,key=access,
    secret=secret,client_kwargs={'endpoint_url': 'https://a3s.fi'})

S3 = boto3.resource(service_name='s3', 
                    aws_access_key_id=access, 
                    aws_secret_access_key=secret, 
                    endpoint_url='https://a3s.fi')


'''





# --- Basic routines ---

def ls_s3(fs, remote_path):
    
    return fs.ls(remote_path)
    

def isfile_s3(fs, remote_path):
    
    return fs.isfile(remote_path)





# --- Reading from S3 ---


def read_mf_s3(fs, remote_paths, **kwargs):
    
    # Open a set of remote netCDF files, leave them open
    opened_files = []
    for file_path in remote_paths:
        opened_files.append(fs.open(file_path, 'rb'))
    
    # Read lazily
    ds = xr.open_mfdataset(opened_files, **kwargs)

    return ds



def read_nc_s3(fs, remote_path, **kwargs):
    
    # Read one file into memory, close after reading
    with fs.open(remote_path, 'rb') as f:
        ds = xr.open_dataset(f, **kwargs)

    return ds


def read_zr_s3(list_of_remote_zarr_files, **kwargs):
    
    # Might not work
    ds = xr.open_mfdataset(list_of_remote_zarr_files, engine='zarr', 
                backend_kwargs=dict(storage_options={'anon': True}))
    
    return ds



# --- Uploading existing data from disk into S3 ---

def upload_data_s3(fs, remote_path, local_path):
    
    fs.put(local_path, remote_path)
    fs.chmod(remote_path, acl='public-read')
    
    pass





# --- Saving special data directly from memory into S3 --- 


def xr_to_zarr_s3(ds, remote_path, **kwargs):
    
    # remote_path = 's3://resiclim/test/2m_temperature/2m_temperature_era5Land_202102.zarr'
    target = fsspec.get_mapper(remote_path,
        anon=False, key=access, secret=secret,
        client_kwargs={'endpoint_url': 'https://a3s.fi'})
    
    task = ds.to_zarr(target, mode='w')
    
    return task




def csv_save_s3(fs, df, remote_path, **kwargs):
    
    with fs.open(remote_path,'w') as f: 
        df.to_csv(f, **kwargs)
    
    fs.chmod(remote_path, acl='public-read')
    
    pass



def fig_save_s3(s3, fig, bucket, remote_path, **kwargs):
    
    img_data = io.BytesIO()
    fig.savefig(img_data, format='png', **kwargs)
    
    s3.Object(bucket,remote_path).put(Body=img_data.getvalue(), ContentType='image/png')
    s3.ObjectAcl(bucket,remote_path).put(ACL='public-read')
    
    pass



