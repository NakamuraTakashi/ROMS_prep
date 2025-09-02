#!/usr/bin/env python3
#---------------------------------------------------------------------------
# Please check the version of Python in which cdsapi is installed 
# by following the command:
#
#> pip show cdsapi
# * check "pythonX.X" written in "Location: /usr/lib/pythonX.X/site-packages"
#
# Executable command (if pythonX.X = python3.12)
#> python3.12 get_era5_atm.py
#---------------------------------------------------------------------------
import os
import cdsapi
#
# Geographical and Grid parameters --------
#
#                 ______ (Rlon,Tlat)
#                |      |
#                |      |
#                |______|
#     (Llon,Blat)                     
#
Llon = 131    # Longitude (degrees) of the bottom-left corner of the grid. 
Rlon = 138    # Longitude (degrees) of the top-right corner of the grid. 
Blat = 4      # Latitude  (degrees) of the bottom-left corner of the grid.
Tlat = 11     # Latitude  (degrees) of the top-right corner of the grid.

#Years=["2017","2018","2019","2020"]
Years=["2022"]

Months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]
#Months = ["01"]

OUTPUT_DIR = 'Data/Palau'

FILE_NAME_prefix='era5_atm'

#---------------------------------------------------------------------------
if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

#-Loop start --------------------------------------------------------------------------

for Year in Years:
    for Month in Months:
    
#        FILE_NAME = OUTPUT_DIR + '/' + FILE_NAME_prefix + '_' + Year + Month + '.nc'
        FILE_NAME = OUTPUT_DIR + '/' + FILE_NAME_prefix + '_' + Year + Month + '.grib'
        
        dataset = "reanalysis-era5-single-levels"
        request = {
            "product_type": ["reanalysis"],
            "variable": [
                "10m_u_component_of_wind",
                "10m_v_component_of_wind",
                "2m_dewpoint_temperature",
                "2m_temperature",
                "mean_sea_level_pressure",
                "total_precipitation",
                "surface_solar_radiation_downwards",
                "surface_thermal_radiation_downwards"
            ],
            "year": Year,
            "month": Month,
            "day": [
                "01", "02", "03",
                "04", "05", "06",
                "07", "08", "09",
                "10", "11", "12",
                "13", "14", "15",
                "16", "17", "18",
                "19", "20", "21",
                "22", "23", "24",
                "25", "26", "27",
                "28", "29", "30", 
                "31"
            ],
            "time": [
                "00:00", "01:00", "02:00",
                "03:00", "04:00", "05:00",
                "06:00", "07:00", "08:00",
                "09:00", "10:00", "11:00",
                "12:00", "13:00", "14:00",
                "15:00", "16:00", "17:00",
                "18:00", "19:00", "20:00",
                "21:00", "22:00", "23:00"
            ],
#            "data_format": "netcdf",
            "data_format": "grib",
            "download_format": "unarchived",
            "area": [Tlat, Llon, Blat, Rlon]
        }
        target = FILE_NAME
        client = cdsapi.Client()
        client.retrieve(dataset, request, target)
        