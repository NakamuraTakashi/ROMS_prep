#!/bin/bash

#gfortran -fbounds-check -fno-align-commons mod_calendar.f90 mod_interpolation.f90 mod_roms_netcdf.f90 set_scoord.f90 naotide.f  pzcon.f potmp.f iniHYCOMpNAO2ROMSv4.F90 -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o iniHYCOMpNAO2ROMS.exe
gfortran -fbounds-check -fno-align-commons mod_calendar.f90 mod_interpolation.f90 mod_roms_netcdf.f90 set_scoord.f90 naotidej.f  pzcon.f potmp.f iniHYCOMpNAO2ROMSv4.F90 -fopenmp -O2 -I/usr/include -L/usr/lib -lnetcdff -o iniHYCOMpNAO2ROMS.exe

export OMP_NUM_THREADS=12

./iniHYCOMpNAO2ROMS.exe
