#!/bin/bash
# If you want to merge HYCOM files, please install CDO (Climate Data Operators).
#  CDO is a command-line tool for processing climate and forecast data.
#  > sudo apt-get install cdo

# The SKIP_SAME_TIME flag means merge files by skipping the overlapping periods.
export SKIP_SAME_TIME=1

#cdo mergetime Palau_HYCOM_extracted_1_20240810.nc \
#              Palau_HYCOM_extracted_1_20240810.nc \
#              Palau_HYCOM_extracted_1_20240810-20250101.nc

#cdo mergetime Palau_HYCOM_extracted_2_20240810.nc \
#              Palau_HYCOM_extracted_2_20240810.nc \
#              Palau_HYCOM_extracted_2_20240810-20250101.nc

#cdo mergetime Palau_HYCOM_extracted_3_20240810.nc \
#              Palau_HYCOM_extracted_3_20240810.nc \
#              Palau_HYCOM_extracted_3_20240810-20250101.nc

#cdo mergetime Palau_HYCOM_extracted_4_20240810.nc \
#              Palau_HYCOM_extracted_4_20240810.nc \
#              Palau_HYCOM_extracted_4_20240810-20250101.nc

#cdo mergetime Palau_HYCOM_extracted_5_20240810.nc \
#              Palau_HYCOM_extracted_5_20240810.nc \
#              Palau_HYCOM_extracted_5_20240810-20250101.nc
#cdo mergetime /mnt/d/HYCOM_DATA/Palau/Palau_HYCOM_extracted_20240208.nc \
#              /mnt/d/HYCOM_DATA/Palau/Palau_HYCOM_extracted_20240211.nc \
#              /mnt/d/HYCOM_DATA/Palau/Palau_HYCOM_extracted_20240215.nc \
#              Palau_HYCOM_extracted_20240208.nc


#cdo merge Palau_HYCOM_extracted_2_20240810.nc \
#          Palau_HYCOM_extracted_3_20240810.nc \
#          Palau_HYCOM_extracted_4_20240810.nc \
#          Palau_HYCOM_extracted_5_20240810.nc \
#          Palau_HYCOM_extracted_2-5_20240810.nc

cdo merge /mnt/d/HYCOM_DATA/Palau/backup/Palau_HYCOM_ESPC_D_V02_2_20250101.nc \
          /mnt/d/HYCOM_DATA/Palau/backup/Palau_HYCOM_ESPC_D_V02_3_20250101.nc \
          /mnt/d/HYCOM_DATA/Palau/backup/Palau_HYCOM_ESPC_D_V02_4_20250101.nc \
          /mnt/d/HYCOM_DATA/Palau/backup/Palau_HYCOM_ESPC_D_V02_5_20250101.nc \
          /mnt/d/HYCOM_DATA/Palau/Palau_HYCOM_ESPC_D_V02_2-5_20250101-20250601.nc
