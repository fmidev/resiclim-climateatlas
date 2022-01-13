#!/bin/sh
#
#
# Download ERA5 Land reanalysis data in CSC Puhti
# Transfer the files immediately to CSC Allas S3 object storage


# Prepare the Anaconda Python environment 
export PATH="/fmi/projappl/project_2002138/miniconda/bin:$PATH"
export LD_LIBRARY_PATH="/fmi/projappl/project_2002138/miniconda/lib"




# Download the variables
#declare -a vars=('2m_dewpoint_temperature' '2m_temperature' 'snow_cover' 'snow_depth' 'snowfall' 'total_precipitation') 
declare -a vars=('2m_temperature' 'total_precipitation') 
for var in "${vars[@]}"
do
   echo $var
   python download_era5_land.py $var &
done



