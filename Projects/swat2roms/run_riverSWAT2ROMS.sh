#!/bin/bash
# ===== Watershed Models ========================================================================
# Please choose one of the following options
#
#export MY_CPP_FLAGS="-DSWAT_PLUS_REV_2019_59_2"
export MY_CPP_FLAGS="-DSWAT_PLUS_REV_2019_59_3"
#export MY_CPP_FLAGS="-DSWAT_PLUS_OKMT"

# Please choose one of the following input file type
#
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DCHANNEL_DAY"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DCHANNEL_SUBDAY"

SRC_DIR=../../src

FCFLAGS="-O2"

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_utility.F90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/riverSWAT2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o riverSWAT2ROMS.exe
rm *.mod

export OMP_NUM_THREADS=12

./riverSWAT2ROMS.exe < test_swat2roms.in

