#!/bin/bash

gfortran mod_utility.f90 mod_interpolation.f90 mod_roms_netcdf.f90 grdGEBCO2ROMS.F90 -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o grdGEBCO2ROMS.exe

export OMP_NUM_THREADS=12

./grdGEBCO2ROMS.exe
