#!/bin/sh

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

#gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcJRA55.F90 -O2 ${INCLUDE} ${LIB} -o frcJRA55.exe
gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcJRA55.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o frcJRA55.exe

export OMP_NUM_THREADS=12
rm *.mod

#./frcJRA55.exe < Berau1.in
./frcJRA55.exe < Berau2.in
