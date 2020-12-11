#!/bin/sh

SRC_DIR=../../src

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/getJMACWM.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o getJMACWM.exe
#gfortran getJMACWMv1.F90 ${INCLUDE} ${LIB} -o getJMACWM.exe

export OMP_NUM_THREADS=12

rm *.mod

./getJMACWM.exe < Shiraho_reef1.in
