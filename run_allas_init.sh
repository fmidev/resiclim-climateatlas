#!/bin/sh -l
# 
#
#
# Commands to initialize the S3 storage system in CSC Allas for the RESICLIM project
#
# Some of these commands work only in the CSC Puhti supercomputer
#
# This file must be run only once in the beginning of the project
#
#
#  S3 Endpoint: a3s.fi
#  DNS-style bucket+hostname:port template for accessing a bucket: %(bucket)s.a3s.fi
#  website_endpoint = http://%(bucket)s.s3-website-%(location)s.amazonaws.com/
#
#  The public address of data (if exists) will follow this format:
#  https://resiclim-era5.a3s.fi/some_folder_path/data.nc




mkdir -p /users/kamarain/resiclim
cd /users/kamarain/resiclim

# Allas environment
module load allas

# Create the configuration file for the s3cmd tool and rename it
allas-conf project_2005404 --mode s3cmd
mv $HOME/.s3cfg $HOME/.s3cfg_resiclim


# Use RESICLIM credentials 
alias s3cmd_resic='s3cmd --config=/users/kamarain/.s3cfg_resiclim'





# Create a bucket for storing the ERA5 files 
s3cmd_resic mb --acl-public   s3://resiclim


# s3cmd_resic put 2m_temperature_era5Land_202001.nc s3://resiclim-era5/2m_temperature/ -P
# https://resiclim.a3s.fi/2m_temperature/2m_temperature_era5Land_202001.nc



# Delete everything
#s3cmd_resic rm  s3://resiclim/2m_temperature/*
#s3cmd_resic rm  s3://resiclim/total_precipitation/*
#s3cmd_resic rb  s3://resiclim


