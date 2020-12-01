#!/bin/bash
# ---- Ocean Models --------------
# Please choose one of the follwing options
#
#export MY_CPP_FLAGS="-DHYCOM"
#export MY_CPP_FLAGS="-DJCOPE"
export MY_CPP_FLAGS="-DROMS"

# ---- For ROMS option --------------
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DWET_DRY"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
#FCFLAGS="-pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.f90  ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/iniOCN2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o iniOCN2ROMS.exe

rm *.mod

export OMP_NUM_THREADS=12

./iniOCN2ROMS.exe < Busuanga1.in
