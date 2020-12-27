#!/bin/sh

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

SRC_DIR=../../src

export MY_CPP_FLAGS="-DANGLE_DIR='${SRC_DIR}'"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DLOCAL_TIME"

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/frcJMAobs2ROMS.F90 ${MY_CPP_FLAGS} -O2 ${INCLUDE} ${LIB} -o frcJMAobs2ROMS.exe

rm *.mod

./frcJMAobs2ROMS.exe < Fukido2.in
