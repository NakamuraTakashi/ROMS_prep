#!/bin/bash
export MY_CPP_FLAGS="-DGOFS_31 -DANALYSIS"
# export MY_CPP_FLAGS="-DGOFS_31 -DREANALYSIS"
# export MY_CPP_FLAGS="-DGOFS_30 -DANALYSIS"
# export MY_CPP_FLAGS="-DGOFS_30 -DREANALYSIS"
#
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSKIP_CHECK_TIME"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/getHYCOM.F90 ${MY_CPP_FLAGS} -O2 -I/usr/include -L/usr/lib -lnetcdff -o getHYCOM.exe
rm *.mod

./getHYCOM.exe < CT_0.08.in
