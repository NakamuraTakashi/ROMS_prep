#!/bin/sh
# SWRAD: Downward short-wave radiation flag
# Since 2017-Dec-05, short wave radiation data has been inclueded in MSM
export MY_CPP_FLAGS="-DSWRAD"
export MY_CPP_FLAGS="-DNETCDF_INPUT"

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcJMAMSM.F90 ${MY_CPP_FLAGS} -fopenmp -O2 ${INCLUDE} ${LIB} -o frcJMAMSM.exe

export OMP_NUM_THREADS=12
rm *.mod

./frcJMAMSM.exe < TokyoBay2.in
