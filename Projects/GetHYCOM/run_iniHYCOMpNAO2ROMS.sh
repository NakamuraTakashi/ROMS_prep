#!/bin/bash

tide_model=naotide
#tide_model=naotidej

SRC_DIR=../../src

echo "tide model: ${tide_model}"

if [[ ${tide_model} == naotide ]]; then
	MY_CPP_FLAGS="-DNAOTIDE";
	MY_CPP_FLAGS="${MY_CPP_FLAGS} -DOMAP_DIR='${SRC_DIR}/omap'"
	gfortran -fbounds-check -fno-align-commons ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.F90 ${SRC_DIR}/naotide.F ${SRC_DIR}/pzcon.f ${SRC_DIR}/potmp.f ${SRC_DIR}/iniHYCOMpNAO2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o iniHYCOMpNAO2ROMS.exe;
elif [[ ${tide_model} == naotidej ]]; then
	MY_CPP_FLAGS="-DNAOTIDEJ";
	MY_CPP_FLAGS="${MY_CPP_FLAGS} -DOMAPJ_DIR='${SRC_DIR}/omapj'";
	gfortran -fbounds-check -fno-align-commons ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/set_scoord.f90 ${SRC_DIR}/set_depth.F90 ${SRC_DIR}/naotidej.F ${SRC_DIR}/pzcon.f ${SRC_DIR}/potmp.f ${SRC_DIR}/iniHYCOMpNAO2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o iniHYCOMpNAO2ROMS.exe;
else
	MY_CPP_FLAGS="-DOMAP_DIR='${SRC_DIR}/omap'";
	gfortran -fbounds-check -fno-align-commons ${SRC_DIR}/mod_calendar.f90 ${SRC_DIR}/mod_interpolation.f90 ${SRC_DIR}/mod_roms_netcdf.F90 ${SRC_DIR}/set_scoord.f90  ${SRC_DIR}/set_depth.F90 ${SRC_DIR}/naotide.F ${SRC_DIR}/pzcon.f ${SRC_DIR}/potmp.f ${SRC_DIR}/iniHYCOMpNAO2ROMS.F90 ${MY_CPP_FLAGS} -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o iniHYCOMpNAO2ROMS.exe;
fi

rm *.mod

export OMP_NUM_THREADS=12

./iniHYCOMpNAO2ROMS.exe < CT_0.08.in
