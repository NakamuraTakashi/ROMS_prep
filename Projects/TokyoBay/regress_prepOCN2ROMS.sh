#!/bin/bash
# ============================================================================================
#  Regression test for the prepOCN2ROMS refactor (branch refactor/prepOCN2ROMS).
#
#  Builds the CURRENT working-tree src/prepOCN2ROMS.F90 AND the committed `master` version,
#  runs both on the same input (regress_TokyoBay2_movejpn.in), and checks the two output
#  netCDF files are byte-identical (cdo diffn + cmp). Use after each refactor step to confirm
#  behaviour is preserved.
#
#  Run from this directory (Projects/TokyoBay) — the .in uses a relative GRID_FILE path that
#  resolves from here.
#
#    ./regress_prepOCN2ROMS.sh [his|bry|ini]      # default: his ;  'all' runs all three
#
#  Requires: gfortran, netCDF-Fortran, cdo, and a git repo (uses `git show master:...`).
#  Machine-specific data paths live in regress_TokyoBay2_movejpn.in (edit there).
# ============================================================================================
set -eu

SRC=../../src
IN=regress_TokyoBay2_movejpn.in
COMMON="$SRC/mod_calendar.f90 $SRC/mod_utility.F90 $SRC/mod_interpolation.f90 \
        $SRC/mod_jcope.F90 $SRC/mod_roms_netcdf.F90 $SRC/set_scoord.f90 \
        $SRC/set_depth.F90 $SRC/pzcon.f $SRC/potmp.f"
FF="-fbounds-check -fno-align-commons -fopenmp -O2"
LINK="-I/usr/include -L/usr/lib -lnetcdff"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

run_one() {            # $1 = his|bry|ini
  local mode=$1 cpp pref
  case $mode in
    his) cpp="-DMOVEJPN_MODEL -DHIS_MODE"; pref=regress_his ;;
    bry) cpp="-DMOVEJPN_MODEL -DBRY_MODE"; pref=regress_bry ;;
    ini) cpp="-DMOVEJPN_MODEL -DINI_MODE"; pref=regress_ini ;;
    *)   echo "unknown mode '$mode' (use his|bry|ini|all)"; return 2 ;;
  esac
  local flags="$cpp -DHYCOM_TIME_DIR=$SRC"
  local work; work=$(mktemp -d)

  echo "==== [$mode] build branch (working tree) + master ===="
  ( rm -f "$work"/*.mod; gfortran $FF $COMMON $SRC/prepOCN2ROMS.F90 $flags $LINK -o "$work/branch.exe" )
  git show master:src/prepOCN2ROMS.F90 > "$work/master_prep.F90"
  ( rm -f "$work"/*.mod; gfortran $FF $COMMON "$work/master_prep.F90" $flags $LINK -o "$work/master.exe" )

  echo "==== [$mode] run master + branch ===="
  rm -f ${pref}_*.nc
  "$work/master.exe" < "$IN" > "$work/master.log" 2>&1 || { echo "[$mode] MASTER run failed"; tail -5 "$work/master.log"; rm -rf "$work"; return 1; }
  mv ${pref}_*.nc "$work/master.nc"
  "$work/branch.exe" < "$IN" > "$work/branch.log" 2>&1 || { echo "[$mode] BRANCH run failed"; tail -5 "$work/branch.log"; rm -rf "$work"; return 1; }
  mv ${pref}_*.nc "$work/branch.nc"

  echo "==== [$mode] compare ===="
  local rc=0
  if cmp -s "$work/master.nc" "$work/branch.nc"; then
    echo "[$mode] PASS: byte-identical to master"
  elif [ -z "$(cdo diffn "$work/master.nc" "$work/branch.nc" 2>/dev/null)" ]; then
    echo "[$mode] PASS: cdo diffn shows 0 differing records (header/byte diffs only)"
  else
    echo "[$mode] FAIL: outputs differ:"; cdo diffn "$work/master.nc" "$work/branch.nc" 2>&1 | tail -5; rc=1
  fi
  rm -rf "$work"
  return $rc
}

modes=${1:-his}
status=0
if [ "$modes" = all ]; then
  for m in his bry ini; do run_one "$m" || status=1; done
else
  run_one "$modes" || status=1
fi
echo "==== regression result: $([ $status -eq 0 ] && echo ALL PASS || echo FAIL) ===="
exit $status
