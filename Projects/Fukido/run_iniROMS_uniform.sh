#!/bin/bash
export MY_CPP_FLAGS="-DWET_DRY"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
#FCFLAGS="-pedantic -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.f90  ${SRC_DIR}/mod_utility.F90 ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/iniROMS_uniform.F90 ${MY_CPP_FLAGS} ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o iniROMS_uniform.exe

rm *.mod

#./iniROMS2ROMS.exe < Fukido20m_parent100m.in
./iniROMS_uniform.exe < Fukido3.in