#!/bin/bash

SRC_DIR=../../src

LIB="-L/usr/lib -L/usr/local/lib -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.f90 ${SRC_DIR}/frcJMAMSM.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o frcJMAMSM.exe
rm *.mod

export OMP_NUM_THREADS=12

./frcJMAMSM.exe < Yaeyama1.in
