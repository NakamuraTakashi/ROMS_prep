#!/bin/sh
#
#rm *.exe

export MY_CPP_FLAGS=""
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DUTM_COORD"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DOUTPUT_SWAN_GRID"

LIB="-L/usr/lib -L/usr/local/lib -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/windROMS2SWAN.F90 ${MY_CPP_FLAGS} -fopenmp -O2 ${INCLUDE} ${LIB} -o windROMS2SWAN.exe

export OMP_NUM_THREADS=11
export HDF5_DISABLE_VERSION_CHECK=1
rm *.mod

./windROMS2SWAN.exe < Palau1.in
