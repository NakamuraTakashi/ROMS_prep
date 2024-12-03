#!/bin/sh
rm frcATM2SWAT.exe
# ---- Atmospheric Models --------------
# Please choose one of the following options
#
#export MY_CPP_FLAGS="-DJMA_MSM"
export MY_CPP_FLAGS="-DDSJRA55"
#export MY_CPP_FLAGS="-DJRA55"
#export MY_CPP_FLAGS="-DERA5"

# ---- For JMA_MSM option --------------
# SWRAD: Downward short-wave radiation flag
# Since 2017-Dec-05, short wave radiation data has been inclueded in MSM
#
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSWRAD"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DNETCDF_INPUT"


LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcATM2SWAT.F90 ${MY_CPP_FLAGS} -fopenmp -O2 ${INCLUDE} ${LIB} -o frcATM2SWAT.exe

export OMP_NUM_THREADS=11
export HDF5_DISABLE_VERSION_CHECK=1
rm *.mod

./frcATM2SWAT.exe < Ishigaki_atm2swat.in
#./frcATM2SWAT.exe < Iriomote_atm2swat.in
