#!/usr/bin/env python3.7
#---------------------------------------------------------------------------
# Please check the version of Python in which cdsapi is installed 
# by following the command:
#
#> pip show cdsapi
# * check "pythonX.X" written in "Location: /usr/lib/pythonX.X/site-packages"
#
# Executable command (if pythonX.X = python3.7)
#> python3.7 get_era5_atm.py
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
Llon = 118    # Longitude (degrees) of the bottom-left corner of the grid. 
Rlon = 127    # Longitude (degrees) of the top-right corner of the grid. 
Blat = 9      # Latitude  (degrees) of the bottom-left corner of the grid.
Tlat = 15     # Latitude  (degrees) of the top-right corner of the grid.

Year='2022'

OUTPUT_DIR = '../../Data/era5/Panay'

FILE_NAME_prefix='era5_atm'

#---------------------------------------------------------------------------
if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

FILE_NAME = OUTPUT_DIR + '/' + FILE_NAME_prefix + '_' + Year + '.nc'

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            '2m_temperature', 'mean_sea_level_pressure', 'surface_solar_radiation_downwards',
            'surface_thermal_radiation_downwards', 'total_precipitation',
        ],
        'year': Year,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            Tlat,  Llon,
            Blat,  Rlon,
        ],
        'format': 'netcdf',
    },
    FILE_NAME )