#!/bin/sh

SRC_DIR=../../src

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

#gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.f90 ${SRC_DIR}/frcJRA55.F90 -O2 ${INCLUDE} ${LIB} -o frcJRA55.exe
gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.f90 ${SRC_DIR}/frcDSJRA55.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o frcDSJRA55.exe
rm *.mod

export OMP_NUM_THREADS=12

./frcDSJRA55.exe < Yaeyama1.in
