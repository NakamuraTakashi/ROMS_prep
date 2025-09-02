#!/bin/bash
#
rm grdROMS_add_ncparams.exe

#cp input/Shizugawa3_grd_v0.2.nc input/Shizugawa3_grd_v0.3b.nc

export MY_CPP_FLAGS=""
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAQUACULTURE1"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAQUACULTURE2"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DUTM_COORD"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DOUTPUT_SWAN_GRID"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-fopenmp -O2"

gfortran ${SRC_DIR}/ll2utm.f ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/grdROMS_add_ncparams.F90 ${MY_CPP_FLAGS} ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o grdROMS_add_ncparams.exe
rm *.mod

export OMP_NUM_THREADS=32

./grdROMS_add_ncparams.exe < Shizugawa3.in
