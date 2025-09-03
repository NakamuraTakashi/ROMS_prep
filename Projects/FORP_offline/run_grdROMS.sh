#!/bin/bash
rm *.exe

export MY_CPP_FLAGS=""
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGEBCO2ROMS"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DFORP_BATH"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DMYBATH2ROMS"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGRID_REFINEMENT"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DUTM_COORD"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBATH_SMOOTHING"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSKIP_ONE_GRID_BAY_REMOVAL"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/utm2ll.f ${SRC_DIR}/mod_utility.F90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/grdROMS.F90 ${MY_CPP_FLAGS} -fopenmp ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o grdROMS.exe
rm *.mod

export OMP_NUM_THREADS=30

./grdROMS.exe < FORP_offline.in
