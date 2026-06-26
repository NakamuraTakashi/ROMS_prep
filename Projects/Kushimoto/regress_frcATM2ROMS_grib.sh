#!/bin/bash
# ============================================================================================
#  Regression test for the frcATM2ROMS refactor (branch refactor/frcATM2ROMS) -- JMA_MSM GRIB.
#
#  Compares the CURRENT working-tree src/frcATM2ROMS.F90 against the committed `master`
#  version (master's JMA_MSM GRIB build compiles AND runs, so no patch is needed) on the same
#  3-hourly MSM GPV input, and checks the two outputs are identical (cmp, with cdo diffn
#  fallback for the julian-round-trip time-ULP difference).
#
#  Run from this directory (Projects/Kushimoto).
#    ./regress_frcATM2ROMS_grib.sh     # JMA_MSM GRIB, 2017-01-01 (8 files x 3 FH = 24 steps)
#
#  Requires: gfortran, netCDF-Fortran, eccodes, cdo, git, and MSM GPV data under ATM_dir.
# ============================================================================================
set -eu

SRC=../../src
IN=regress_Kushimoto_msmgrib.in
CPP="-DJMA_MSM"
COMMON="$SRC/mod_calendar.f90 $SRC/mod_interpolation.f90 $SRC/mod_roms_netcdf.F90 $SRC/mod_infile.F90"
FF="-fbounds-check -fno-align-commons -fopenmp -O2"
LINK="-I/usr/include -I/usr/local/include -L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

DATETAG=20170101    # start date 2017-01-01 -> _20170101_{1,2}.nc

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

sed "s#REGRESS_FRC_PREFIX_PLACEHOLDER#$work/out#g" "$IN" > "$work/run.in"

echo "==== build branch (working tree) + master ===="
rm -f "$work"/*.mod
gfortran $FF $COMMON $SRC/frcATM2ROMS.F90 $CPP $LINK -o "$work/branch.exe"
git show master:src/frcATM2ROMS.F90 > "$work/master_frc.F90"
rm -f "$work"/*.mod
gfortran $FF $COMMON "$work/master_frc.F90" $CPP $LINK -o "$work/master.exe"

run_exe() { local exe=$1 tag=$2
  ( cd "$(dirname "$0")" && "$exe" < "$work/run.in" > "$work/$tag.log" 2>&1 ) \
    || { echo "[$tag] run failed"; tail -8 "$work/$tag.log"; exit 1; }
}

echo "==== run master ===="
run_exe "$work/master.exe" master
mv "$work"/out_${DATETAG}_1.nc "$work/master_1.nc"
mv "$work"/out_${DATETAG}_2.nc "$work/master_2.nc"

echo "==== run branch ===="
run_exe "$work/branch.exe" branch
mv "$work"/out_${DATETAG}_1.nc "$work/branch_1.nc"
mv "$work"/out_${DATETAG}_2.nc "$work/branch_2.nc"

echo "==== compare ===="
rc=0
for part in 1 2; do
  m="$work/master_${part}.nc" ; b="$work/branch_${part}.nc"
  if cmp -s "$m" "$b"; then
    echo "[_${part}.nc] PASS: byte-identical to master"
  elif [ -z "$(cdo diffn "$m" "$b" 2>/dev/null)" ]; then
    echo "[_${part}.nc] PASS: cdo diffn 0 differing records (header/time-ULP diffs only)"
  else
    echo "[_${part}.nc] FAIL: outputs differ:"; cdo diffn "$m" "$b" 2>&1 | tail -8; rc=1
  fi
done
exit $rc
