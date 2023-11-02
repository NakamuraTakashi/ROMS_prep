#!/bin/bash

export HDF5_DISABLE_VERSION_CHECK=1

SCRIP_DIR=../../SCRIP_COAWST

mpirun -use-hwthread-cpus -np 10 ${SCRIP_DIR}/scrip_coawst.exe scrip_coawst_shiraho_reef.in 
#