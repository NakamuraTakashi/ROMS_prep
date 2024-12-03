#!/bin/sh
#
rm grdROMS2SWAN.exe

export MY_CPP_FLAGS=""
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DUTM_COORD"

LIB="-L/usr/lib -L/usr/local/lib -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/grdROMS2SWAN.F90 ${MY_CPP_FLAGS} -fopenmp -O2 ${INCLUDE} ${LIB} -o grdROMS2SWAN.exe

export OMP_NUM_THREADS=7
export HDF5_DISABLE_VERSION_CHECK=1
rm *.mod

#./grdROMS2SWAN.exe < Shizugawa2.in
./grdROMS2SWAN.exe < Shizugawa3.in
