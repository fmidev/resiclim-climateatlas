#!/bin/sh
#
#
# Download ERA5 Land reanalysis data in CSC Puhti
# Transfer the files immediately to CSC Allas S3 object storage




code_dir='/users/kamarain/resiclim-climateatlas'
data_dir='/scratch/project_2005404'



# Prepare the Anaconda Python environment 
export PATH="/projappl/project_2005404/miniconda/bin:$PATH"
export LD_LIBRARY_PATH="/projappl/project_2005404/miniconda/lib"

#export PATH="/fmi/projappl/project_2002138/miniconda/bin:$PATH"
#export LD_LIBRARY_PATH="/fmi/projappl/project_2002138/miniconda/lib"



mkdir -p $data_dir/temporal
cd $data_dir/temporal


# Download the variables
#declare -a vars=('2m_dewpoint_temperature' '2m_temperature' 'skin_temperature' 'snow_cover' 'snow_depth' 'snowfall' 'total_precipitation') 
declare -a vars=('2m_temperature' 'total_precipitation') 
for var in "${vars[@]}"
do
   echo $var
   python $code_dir/download_era5_land_zr.py $var &
done



