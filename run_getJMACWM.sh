#!/bin/sh

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

gfortran getJMACWMv1.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o getJMACWM.exe
#gfortran getJMACWMv1.F90 ${INCLUDE} ${LIB} -o getJMACWM.exe

export OMP_NUM_THREADS=12

./getJMACWM.exe
