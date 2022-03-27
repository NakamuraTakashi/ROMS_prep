#!/bin/bash

export MY_CPP_FLAGS="-DNAOTIDEJ"
#export MY_CPP_FLAGS="-DNAOTIDE"

# ========================================================================================

SRC_DIR=../../src

MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSRC_DIR='${SRC_DIR}'"

gfortran -fbounds-check -fno-align-commons ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_utility.F90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_tide.F90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcTIDE2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o frcTIDE2ROMS.exe;

rm *.mod

export OMP_NUM_THREADS=12

./frcTIDE2ROMS.exe < TokyoBay1.in
