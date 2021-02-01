#!/bin/sh
# ---- Atmospheric Models --------------
# Please choose one of the following options
#
export MY_CPP_FLAGS="-DJMA_MSM"
#export MY_CPP_FLAGS="-DDSJRA55"
#export MY_CPP_FLAGS="-DJRA55"

# ---- For JMA_MSM option --------------
# SWRAD: Downward short-wave radiation flag
# Since 2017-Dec-05, short wave radiation data has been inclueded in MSM
#
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSWRAD"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DNETCDF_INPUT"

# ---- For JRA55 or DSJRA55 option --------------
# BULK_FLUX: Downward/Upward short-wave/long-wave radiation flag
# 
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBULK_FLUX"

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcATM2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp -O2 ${INCLUDE} ${LIB} -o frcATM2ROMS_nc.exe

export OMP_NUM_THREADS=12
rm *.mod

#./frcATM2ROMS_2017.exe < Yaeyama1.in
#./frcATM2ROMS.exe < Yaeyama2.in
#./frcATM2ROMS.exe < Yaeyama3.in
OMP_NUM_THREADS=12 ./frcATM2ROMS_nc.exe < Yaeyama3.in
