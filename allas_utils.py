#!/usr/bin/env python



import io, os
import s3fs, boto3
import xarray as xr
import pandas as pd




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








# --- Basic routines ---

def ls_s3(remote_path):
    
    return FS.ls(remote_path)
    

def isfile_s3(remote_path):
    
    return FS.isfile(remote_path)





# --- Reading from S3 ---

def read_nc_s3(remote_path, **kwargs):
    
    with FS.open(remote_path, 'rb') as f:
        ds = xr.open_dataset(f, **kwargs)

    return ds





# --- Uploading existing data from disk into S3 ---

def upload_data_s3(remote_path, local_path):
    
    FS.put(local_path, remote_path)
    FS.chmod(remote_path, acl='public-read')
    
    pass





# --- Saving special data directly from memory into S3 --- 

def csv_save_s3(df, remote_path, **kwargs):
    
    with FS.open(remote_path,'w') as f: 
        df.to_csv(f, **kwargs)
    
    FS.chmod(remote_path, acl='public-read')
    
    pass



def fig_save_s3(fig, bucket, remote_path, **kwargs):
    
    img_data = io.BytesIO()
    fig.savefig(img_data, format='png', **kwargs)
    
    S3.Object(bucket,remote_path).put(Body=img_data.getvalue(), ContentType='image/png')
    S3.ObjectAcl(bucket,remote_path).put(ACL='public-read')
    
    pass



