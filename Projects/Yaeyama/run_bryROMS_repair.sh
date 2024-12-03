#!/bin/bash
rm bryROMS_repair.exe

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran ${SRC_DIR}/mod_utility.F90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/bryROMS_repair.F90 ${MY_CPP_FLAGS} ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o bryROMS_repair.exe
rm *.mod

./bryROMS_repair.exe < Y1_bryROMS_repair.in

