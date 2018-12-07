#!/bin/bash

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.f90 ${SRC_DIR}/grdGEBCO2ROMS.F90 -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o grdGEBCO2ROMS.exe
rm *.mod

export OMP_NUM_THREADS=12

./grdGEBCO2ROMS.exe < CT_0.08.in
