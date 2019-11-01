#!/bin/bash
export MY_CPP_FLAGS="-DWET_DRY"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
#FCFLAGS="-pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.f90  ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/bryROMS2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o bryROMS2ROMS.exe

rm *.mod

export OMP_NUM_THREADS=1

#./bryROMS2ROMS.exe < Berau1.in
./bryROMS2ROMS.exe < Berau2.in