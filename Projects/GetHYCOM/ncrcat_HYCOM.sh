#!/bin/bash
# If you want to merge HYCOM files, please install CDO (Climate Data Operators).
#  CDO is a command-line tool for processing climate and forecast data.
#  > sudo apt-get install cdo

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


ncecat Palau_HYCOM_extracted_2_20240810.nc \
       Palau_HYCOM_extracted_3_20240810.nc \
       Palau_HYCOM_extracted_4_20240810.nc \
       Palau_HYCOM_extracted_5_20240810.nc \
       Palau_HYCOM_extracted_2-5_20240810.nc

#cdo merge Palau_HYCOM_extracted_2_20240810-20250101.nc \
#          Palau_HYCOM_extracted_3_20240810-20250101.nc \
#          Palau_HYCOM_extracted_4_20240810-20250101.nc \
#          Palau_HYCOM_extracted_5_20240810-20250101.nc \
#          Palau_HYCOM_extracted_2-5_20240810-20250101.nc
