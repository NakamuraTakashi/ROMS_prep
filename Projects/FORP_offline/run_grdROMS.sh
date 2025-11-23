#!/bin/bash
rm *.exe

export MY_CPP_FLAGS=""
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGEBCO2ROMS"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DMOVEJPN_BATH"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DMYBATH2ROMS"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGRID_REFINEMENT"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DUTM_COORD"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBATH_SMOOTHING"
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSKIP_ONE_GRID_BAY_REMOVAL"

SRC_DIR=../../src
#FCFLAGS="-Wall -pedantic -std=f95 -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace"
FCFLAGS="-O2"

gfortran \
  ${SRC_DIR}/utm2ll.f \
  ${SRC_DIR}/mod_utility.F90 \
  ${SRC_DIR}/mod_interpolation.f90 \
  ${SRC_DIR}/mod_roms_netcdf.F90 \
  ${SRC_DIR}/mod_movejpn.F90 \
  ${SRC_DIR}/grdROMS.F90 \
  ${MY_CPP_FLAGS} -fopenmp \
  ${FCFLAGS} -I/usr/include -L/usr/lib -lnetcdff -o grdROMS.exe
rm *.mod

export OMP_NUM_THREADS=30

./grdROMS.exe < FORP_offline.in

# memo: command for trim the specific area
##ncks -d xi_rho,128,779 -d eta_rho,104,791 -d xi_u,128,778 -d eta_u,104,791 -d xi_v,128,779 -d eta_v,104,790 -d xi_psi,128,778 -d eta_psi,104,790 ../../../COAWST_DATA/FORP_offline/Grid/forp-jpn-v4_grd_v1.1.nc -o forp-jpn-south_grd_v1.1.nc
#ncks -d xi_rho,129,866 -d eta_rho,80,856 -d xi_u,129,865 -d eta_u,80,856 -d xi_v,129,866 -d eta_v,80,855 -d xi_psi,129,865 -d eta_psi,80,855 forp-jpn-v4_grd_v2.0.nc -o forp-jpn-south_grd_v2.0.nc
