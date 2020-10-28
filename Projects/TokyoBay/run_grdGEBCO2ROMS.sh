#!/bin/bash
export MY_CPP_FLAGS="-DGRID_REFINEMENT"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/grdGEBCO2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o grdGEBCO2ROMS.exe
rm *.mod

export OMP_NUM_THREADS=6

#./grdGEBCO2ROMS.exe < TokyoBay1.in
#./grdGEBCO2ROMS.exe < TokyoBay2.in
./grdGEBCO2ROMS.exe < TokyoBay3.in
