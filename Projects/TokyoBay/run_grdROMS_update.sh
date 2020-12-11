#!/bin/bash
export MY_CPP_FLAGS=""
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGRID_FINE2COARSE"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBATH_SMOOTHING"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/grdROMS_update.F90 ${MY_CPP_FLAGS} ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o grdROMS_update.exe
rm *.mod

#./grdROMS_update.exe < TokyoBay2.in
./grdROMS_update.exe < TokyoBay2.in
#./grdROMS_update.exe < TokyoBay3.in
