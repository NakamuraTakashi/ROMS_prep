#!/bin/sh
## Run in current working directory
#$ -cwd
## Resource type F: qty 2
#$ -l f_node=1
## maximum run time
#$ -l h_rt=5:00:00
#$ -N bryOCN2ROMS
## Initialize module command
. /etc/profile.d/modules.sh

module load cuda openmpi
module load hdf5-parallel
module load netcdf-parallel

export OMP_NUM_THREADS=28

./bryOCN2ROMS < Yaeyama3_T3.in
#./iniOCN2ROMS < Yaeyama3_T3.in
