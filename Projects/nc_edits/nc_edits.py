#!/usr/bin/env python3
# Need 
# > pip install netcdf4

import netCDF4

NC_FILE = '../Shizugawa/input/Shizugawa3_grd_v0.3b.nc'
nc = netCDF4.Dataset(NC_FILE, 'r')
#dim = len(nc.dimensions['dimname'])
var = nc.variables['aquaculture_01'][:]

nc.close()
        