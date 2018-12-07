#!/bin/bash

LIB="-L/usr/lib -L/usr/local/lib -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

gfortran mod_calendar.f90 mod_utility.f90 mod_interpolation.f90 mod_roms_netcdf.f90 frcJMAMSM.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o frcJMAMSM.exe

export OMP_NUM_THREADS=12

./frcJMAMSM.exe < Yaeyama.in
