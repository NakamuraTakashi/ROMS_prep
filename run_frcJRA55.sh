#!/bin/sh

LIB="-L/usr/lib -lnetcdff -L${HOME}/grib2/lib -lwgrib2 -lgfortran -lz -lm"
MOD="-I/usr/include -I${HOME}/grib2/lib"

gfortran mod_calendar.f90 mod_interpolation.f90 frcJRA55.F90 -O2 ${MOD} ${LIB} -o frcJRA55.exe

./frcJRA55.exe
