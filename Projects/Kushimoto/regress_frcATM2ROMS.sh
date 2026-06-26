#!/bin/bash
# ============================================================================================
#  Regression test for the frcATM2ROMS refactor (branch refactor/frcATM2ROMS) -- JMA_MSM NetCDF.
#
#  Compares the CURRENT working-tree src/frcATM2ROMS.F90 against a *reference* built from the
#  committed `master` version. NOTE: master's JMA_MSM+NETCDF_INPUT is doubly broken:
#    (1) compile: amask/amaskv are declared only for ERA5/FORP;
#    (2) runtime: the netcdf branch writes IN_FILE2(1:2) which is never allocated for netcdf.
#  The reference applies the MINIMAL no-op fixes (declare amask/amaskv for JMA_MSM, and drop
#  the two dead write-only IN_FILE2 lines in the netcdf branch) so it runs the true legacy
#  algorithm. Both fixes leave the output unchanged.
#
#  Run from this directory (Projects/Kushimoto).
#    ./regress_frcATM2ROMS.sh          # JMA_MSM NetCDF, 2010-06-01..06-03 (3 daily files)
#
#  Requires: gfortran, netCDF-Fortran, eccodes, cdo, git, and the MSM data under ATM_dir.
# ============================================================================================
set -eu

SRC=../../src
IN=regress_Kushimoto_msmnc.in
CPP="-DJMA_MSM -DNETCDF_INPUT"
COMMON="$SRC/mod_calendar.f90 $SRC/mod_interpolation.f90 $SRC/mod_roms_netcdf.F90 $SRC/mod_infile.F90"
FF="-fbounds-check -fno-align-commons -fopenmp -O2"
LINK="-I/usr/include -I/usr/local/include -L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

DATETAG=20100601    # start date 2010-06-01 -> _20100601_{1,2}.nc

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

sed "s#REGRESS_FRC_PREFIX_PLACEHOLDER#$work/out#g" "$IN" > "$work/run.in"

echo "==== build branch (working tree) ===="
rm -f "$work"/*.mod
gfortran $FF $COMMON $SRC/frcATM2ROMS.F90 $CPP $LINK -o "$work/branch.exe"

echo "==== build reference (master + minimal amask-declaration fix) ===="
git show master:src/frcATM2ROMS.F90 > "$work/master_frc.F90"
# Minimal compile-only fix: declare the ERA5/FORP-only INFILE vars (incl. amask/amaskv)
# for JMA_MSM as well. This is the single declaration-block guard (the one right before
# `TYPE T_NC`); no executable line changes.
perl -0pi -e 's/#if defined ERA5 \|\| defined FORP_ATM\n  TYPE T_NC/#if defined ERA5 || defined FORP_ATM || defined JMA_MSM\n  TYPE T_NC/' "$work/master_frc.F90"
# Drop the two dead write-only IN_FILE2 lines in the netcdf branch (runtime-crash fix, no-op).
perl -0pi -e 's/    ihours = ihours \+ 24  !!! Files exist daily interval\n\n    IN_FILE2\(1\) = IN_FILE\(1\)  \n    IN_FILE2\(2\) = IN_FILE\(2\)  \n/    ihours = ihours + 24  !!! Files exist daily interval\n\n/' "$work/master_frc.F90"
rm -f "$work"/*.mod
gfortran $FF $COMMON "$work/master_frc.F90" $CPP $LINK -o "$work/master.exe"

run_exe() { local exe=$1 tag=$2
  ( cd "$(dirname "$0")" && "$exe" < "$work/run.in" > "$work/$tag.log" 2>&1 ) \
    || { echo "[$tag] run failed"; tail -8 "$work/$tag.log"; exit 1; }
}

echo "==== run reference ===="
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
    echo "[_${part}.nc] PASS: byte-identical to reference"
  elif [ -z "$(cdo diffn "$m" "$b" 2>/dev/null)" ]; then
    echo "[_${part}.nc] PASS: cdo diffn 0 differing records (header/time-ULP diffs only)"
  else
    echo "[_${part}.nc] FAIL: outputs differ:"; cdo diffn "$m" "$b" 2>&1 | tail -8; rc=1
  fi
done
exit $rc
