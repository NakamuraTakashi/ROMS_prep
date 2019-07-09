#!/bin/bash

SRC_DIR=../../src

#gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.F90 ${SRC_DIR}/pzcon.f ${SRC_DIR}/potmp.f ${SRC_DIR}/bryHYCOM2ROMS.F90 -O2 -I/usr/include -L/usr/lib -lnetcdff -o bryHYCOM2ROMS.exe
gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.F90 ${SRC_DIR}/pzcon.f ${SRC_DIR}/potmp.f ${SRC_DIR}/bryHYCOM2ROMS.F90 -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o bryHYCOM2ROMS.exe
rm *.mod

export OMP_NUM_THREADS=12

./bryHYCOM2ROMS.exe < CT_0.08.in
