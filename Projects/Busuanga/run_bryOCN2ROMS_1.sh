#!/bin/bash
# ===== Ocean Models ========================================================================
# Please choose one of the following options
#
export MY_CPP_FLAGS="-DHYCOM_MODEL"
#export MY_CPP_FLAGS="-DJCOPE_MODEL"
#export MY_CPP_FLAGS="-DROMS_MODEL"

# ===== HYCOM option ====================================================================

# Please choose one of the following options
# -- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2018-12-04 to present)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_31 -DANALYSIS_Y"

# -- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2014-07-01 to 2020-02-18)
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_31 -DANALYSIS"

# -- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Renalysis (since 1994-01-01 to 2015-12-31)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_31 -DREANALYSIS"

# --  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Analysis (since 2008-09-19 to 2018-11-20)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_30 -DANALYSIS"

# --  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Reanalysis (since 1992-10-02 to 2012-12-31)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_30 -DREANALYSIS"

#----------------------------------------------------------------------------------------
# Please activate if you want to skip time checking 
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSKIP_CHECK_TIME"

# Fast read option for HYCOM 4D data (u, v, temp, salt).
# *If failure frequently occurs, please deactivate this option. 
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DFAST_READ"

# ===== ROMS option ======================================================================
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DWET_DRY"

# ===== JCOPE option =====================================================================
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DJCOPE_T"


echo "${MY_CPP_FLAGS}"

SRC_DIR=../../src

export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DHYCOM_TIME_DIR='${SRC_DIR}'"

gfortran -fbounds-check -fno-align-commons ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_jcope.F90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/set_scoord.f90  ${SRC_DIR}/set_depth.F90 ${SRC_DIR}/pzcon.f ${SRC_DIR}/potmp.f ${SRC_DIR}/bryOCN2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o bryOCN2ROMS.exe
rm *.mod

export OMP_NUM_THREADS=12

./bryOCN2ROMS.exe < Berau1.in
