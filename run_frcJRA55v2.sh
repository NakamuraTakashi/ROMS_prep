#!/bin/sh

#LIB="-L/usr/lib -lnetcdff -L${HOME}/grib2/lib -lwgrib2 -lgfortran -lz -lm"
#MOD="-I/usr/include -I${HOME}/grib2/lib"
LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

#gfortran mod_calendar.f90 mod_interpolation.f90 mod_roms_netcdf.f90 frcJRA55v2.F90 -O2 ${INCLUDE} ${LIB} -o frcJRA55.exe
gfortran mod_calendar.f90 mod_interpolation.f90 mod_roms_netcdf.f90 frcJRA55v2.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o frcJRA55.exe

export OMP_NUM_THREADS=12

./frcJRA55.exe
