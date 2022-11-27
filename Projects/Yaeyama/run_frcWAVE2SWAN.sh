#!/bin/sh
# ---- Wave Models --------------
# Please choose one of the following options
#
export MY_CPP_FLAGS="-DJMA_CWM"
#export MY_CPP_FLAGS="-DWW3"

SRC_DIR=../../src

LIB="-L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
INCLUDE="-I/usr/include -I/usr/local/include"

gfortran ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_utility.f90 ${SRC_DIR}/frcWAVE2SWAN.F90 -fopenmp -O2 ${INCLUDE} ${LIB} -o frcWAVE2SWAN.exe

export OMP_NUM_THREADS=12

rm *.mod

./frcWAVE2SWAN.exe < Yaeyama1.in
