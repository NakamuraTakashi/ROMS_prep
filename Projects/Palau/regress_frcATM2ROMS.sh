#!/bin/bash
# ============================================================================================
#  Regression test for the frcATM2ROMS refactor (branch refactor/frcATM2ROMS).
#
#  Builds the CURRENT working-tree src/frcATM2ROMS.F90 AND the committed `master` version,
#  runs both on the same ERA5 (GRIB) input (regress_Palau2_era5.in), and checks the two
#  output netCDF files are identical (cmp, with cdo diffn fallback). Use after each refactor
#  step to confirm behaviour is preserved.
#
#  Run from this directory (Projects/Palau) — the .in uses a relative GRID_FILE and ERA5
#  data path that resolve from here.
#
#    ./regress_frcATM2ROMS.sh                 # ERA5 GRIB, 2023-12-28 .. 2024-01-03 (file switch)
#
#  Requires: gfortran, netCDF-Fortran, eccodes, cdo, and a git repo (uses `git show master:`).
# ============================================================================================
set -eu

SRC=../../src
IN=regress_Palau2_era5.in
CPP="-DERA5"
COMMON="$SRC/mod_calendar.f90 $SRC/mod_interpolation.f90 $SRC/mod_roms_netcdf.F90 $SRC/mod_infile.F90"
FF="-fbounds-check -fno-align-commons -fopenmp -O2"
LINK="-I/usr/include -I/usr/local/include -L/usr/lib -L/usr/local/lib -leccodes_f90 -lnetcdff"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

# Output file prefix (start date 2023-12-28 -> _20231228_{1,2}.nc)
DATETAG=20231228

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# Materialise the .in with a concrete (scratchpad) output prefix.
sed "s#REGRESS_FRC_PREFIX_PLACEHOLDER#$work/out#g" "$IN" > "$work/run.in"

echo "==== build branch (working tree) + master ===="
rm -f "$work"/*.mod
gfortran $FF $COMMON $SRC/frcATM2ROMS.F90 $CPP $LINK -o "$work/branch.exe"
git show master:src/frcATM2ROMS.F90 > "$work/master_frc.F90"
rm -f "$work"/*.mod
gfortran $FF $COMMON "$work/master_frc.F90" $CPP $LINK -o "$work/master.exe"

run_exe() {            # $1 = exe  $2 = tag
  local exe=$1 tag=$2
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
    echo "[_${part}.nc] PASS: cdo diffn 0 differing records (header/byte diffs only)"
  else
    echo "[_${part}.nc] FAIL: outputs differ:"; cdo diffn "$m" "$b" 2>&1 | tail -8; rc=1
  fi
done
exit $rc
