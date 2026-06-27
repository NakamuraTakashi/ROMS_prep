#!/bin/sh
#rm *.exe
# ===== Input namelist file =================================================================
INPUT=TokyoBay2.in

# ---- Atmospheric Models --------------
# Please choose one of the following options
#
export MY_CPP_FLAGS="-DJMA_MSM"
#export MY_CPP_FLAGS="-DDSJRA55"
#export MY_CPP_FLAGS="-DJMA_LSM"
#export MY_CPP_FLAGS="-DJRA55"
#export MY_CPP_FLAGS="-DERA5"
#export MY_CPP_FLAGS="-DFORP_ATM"

# ---- For JMA_MSM option --------------
# SWRAD: Downward short-wave radiation flag
# Since 2017-Dec-05, short wave radiation data has been inclueded in MSM
#
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSWRAD"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DNETCDF_INPUT"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DJMA_MSM_CLOUD_ONLY"

# ---- For JRA55 or DSJRA55 option --------------
# BULK_FLUX: Downward/Upward short-wave/long-wave radiation flag
# 
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBULK_FLUX"

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I${PWD} -I/usr/include -I/usr/local/include"
SRC_DIR=../../src
FFLAGS="-fopenmp -O3"
#FFLAGS="-fbounds-check -ffree-form -O0 -g -fcheck=array-temps,bounds,do,mem,pointer,recursion"

gfortran ${FFLAGS} \
  ${SRC_DIR}/mod_calendar.f90 \
  ${SRC_DIR}/mod_infile.F90 \
  ${SRC_DIR}/mod_interpolation.f90 \
  ${SRC_DIR}/mod_roms_netcdf.F90 \
  ${SRC_DIR}/frcATM2ROMS.F90 \
  ${MY_CPP_FLAGS} ${INCLUDE} ${LIB} -o frcATM2ROMS.exe

export OMP_NUM_THREADS=12
rm *.mod

./frcATM2ROMS.exe < ${INPUT}
