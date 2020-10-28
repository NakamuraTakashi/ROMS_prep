#!/bin/bash
#export MY_CPP_FLAGS="-DGOFS_31 -DANALYSIS_Y"  # since 2018-12-04 to present, 3 hourly
 export MY_CPP_FLAGS="-DGOFS_31 -DANALYSIS"    # since 2014-07-01 to 2020-02-18, 3 hourly
#export MY_CPP_FLAGS="-DGOFS_31 -DREANALYSIS"  # since 1994-01-01 to 2015-12-31, 3 hourly
#export MY_CPP_FLAGS="-DGOFS_30 -DANALYSIS"    # since 2008-09-19 to 2018-11-20, daily
#export MY_CPP_FLAGS="-DGOFS_30 -DREANALYSIS"  # since 1992-10-02 to 2012-12-31, daily
#
 export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSKIP_CHECK_TIME"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/getHYCOM.F90 ${MY_CPP_FLAGS} -O2 -I/usr/include -L/usr/lib -lnetcdff -o getHYCOM.exe
rm *.mod

./getHYCOM.exe < TokyoBay1.in
