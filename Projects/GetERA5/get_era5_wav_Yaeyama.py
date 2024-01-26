#!/usr/bin/env python3.9
#---------------------------------------------------------------------------
# Please check the version of Python in which cdsapi is installed 
# by following the command:
#
#> pip show cdsapi
# * check "pythonX.X" written in "Location: /usr/lib/pythonX.X/site-packages"
#
# Executable command (if pythonX.X = python3.9)
#> python3.9 get_era5_wav.py
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

Llon = 120    # Longitude (degrees) of the bottom-left corner of the grid. 
Rlon = 128    # Longitude (degrees) of the top-right corner of the grid. 
Blat = 21     # Latitude  (degrees) of the bottom-left corner of the grid.
Tlat = 27     # Latitude  (degrees) of the top-right corner of the grid.

Year='2008'

OUTPUT_DIR = 'D:/ERA5/Yaeyama'

FILE_NAME_prefix='ERA5_wav_Yaeyama'

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
            'mean_wave_direction', 
            'mean_wave_period', 
            'significant_height_of_combined_wind_waves_and_swell',
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