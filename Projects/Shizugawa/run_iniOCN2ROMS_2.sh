#!/bin/bash
# ============================================================================================
#  prepOCN2ROMS : unified generator for ROMS initial / boundary / his input files.
#  Select the output with gen_mode below (replaces run_iniOCN2ROMS_*.sh / run_bryOCN2ROMS_*.sh).
# ============================================================================================

# ===== Generation mode =====================================================================
# Choose ONE: ini = initial conditions (single time-step, + tide option)
#             bry = boundary conditions (time series)
#             his = ROMS-his-like file  (time series, whole domain)
#gen_mode=bry
#gen_mode=his
gen_mode=ini

# ===== Input namelist file =================================================================
INPUT=Shizugawa2.in

# ===== Ocean Models ========================================================================
# Please choose one of the following options
#
#export MY_CPP_FLAGS="-DHYCOM_MODEL"
#export MY_CPP_FLAGS="-DJCOPE_MODEL"
export MY_CPP_FLAGS="-DROMS_MODEL"
#export MY_CPP_FLAGS="-DFORP_MODEL"
#export MY_CPP_FLAGS="-DMOVEJPN_MODEL"
#export MY_CPP_FLAGS="-DFORA_MODEL"

# ===== HYCOM option ========================================================================
# Please choose one of the following options
# -- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2018-12-04 to present)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_31 -DANALYSIS_Y"
# -- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Analysis (since 2014-07-01 to 2020-02-18)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_31 -DANALYSIS"
# -- GOFS 3.1: 41-layer HYCOM + NCODA Global 1/12 deg Renalysis (since 1994-01-01 to 2015-12-31)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_31 -DREANALYSIS"
# --  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Analysis (since 2008-09-19 to 2018-11-20)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_30 -DANALYSIS"
# --  GOFS 3.0: HYCOM + NCODA Global 1/12 deg Reanalysis (since 1992-10-02 to 2012-12-31)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DGOFS_30 -DREANALYSIS"
# --  Local HYCOM extracted data by getHYCOM code
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DHYCOM_LOCAL"

#--------------------------------------------------------------------------------------------
# Please activate if you want to skip time checking (recommended for ini)
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DSKIP_CHECK_TIME"

# Fast read option for HYCOM 4D data (u, v, temp, salt).
# *If failure frequently occurs, please deactivate this option.
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DFAST_READ"

# ===== ROMS option =========================================================================
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DWET_DRY"

# ===== JCOPE option ========================================================================
#export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DJCOPE_T"

# ===== Tide model (ini only) ===============================================================
# Adds NAOTIDE/NAOTIDEJ tide to zeta. Leave as 'none' for bry/his (ignored there).
#tide_model=naotide
#tide_model=naotidej
tide_model=none

# ===========================================================================================
SRC_DIR=../../src

# ----- Map gen_mode -> output-mode CPP flag + exe name -----
case ${gen_mode} in
  ini) MODE_FLAG="-DINI_MODE"; EXE=iniOCN2ROMS_3.exe ;;
  bry) MODE_FLAG="-DBRY_MODE"; EXE=bryOCN2ROMS_3.exe ;;
  his) MODE_FLAG="-DHIS_MODE"; EXE=hisOCN2ROMS_3.exe ;;
  *)   echo "ERROR: gen_mode must be one of ini / bry / his"; exit 1 ;;
esac
export MY_CPP_FLAGS="${MY_CPP_FLAGS} ${MODE_FLAG}"
export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DHYCOM_TIME_DIR='${SRC_DIR}'"

# ----- Tide sources / flags (only used for ini) -----
TIDE_SRC=""
if [[ ${gen_mode} == ini ]]; then
  if [[ ${tide_model} == naotide ]]; then
    MY_CPP_FLAGS="${MY_CPP_FLAGS} -DNAOTIDE -DOMAP_DIR='${SRC_DIR}/omap'"
    TIDE_SRC="${SRC_DIR}/naotide.F"
  elif [[ ${tide_model} == naotidej ]]; then
    MY_CPP_FLAGS="${MY_CPP_FLAGS} -DNAOTIDEJ -DOMAPJ_DIR='${SRC_DIR}/omapj'"
    TIDE_SRC="${SRC_DIR}/naotidej.F"
  fi
fi

echo "gen_mode   : ${gen_mode}  (${MODE_FLAG})"
echo "input      : ${INPUT}"
echo "tide_model : ${tide_model}"
echo "CPP flags  : ${MY_CPP_FLAGS}"

# ----- Compile -----
gfortran \
  ${SRC_DIR}/mod_calendar.f90 \
  ${SRC_DIR}/mod_utility.F90 \
  ${SRC_DIR}/mod_infile.F90 \
  ${SRC_DIR}/mod_interpolation.f90 \
  ${SRC_DIR}/mod_jcope.F90 \
  ${SRC_DIR}/mod_roms_netcdf.F90 \
  ${SRC_DIR}/set_scoord.f90 \
  ${SRC_DIR}/set_depth.F90 \
  ${TIDE_SRC} \
  ${SRC_DIR}/pzcon.f \
  ${SRC_DIR}/potmp.f \
  ${SRC_DIR}/prepOCN2ROMS.F90 \
  ${MY_CPP_FLAGS} -fopenmp -O3 -I/usr/include -L/usr/lib -lnetcdff -o ${EXE} || { echo "Compile failed"; exit 1; }

rm -f *.mod

# ----- Run -----
export OMP_NUM_THREADS=12
./${EXE} < ${INPUT}

# NOTE on ini: the single initial time is taken from the &sdate block of the input
#   (the shared engine reuses the his time setup; &ini supplies only INI_prefix).
#   For an initial file at a specific date, set &sdate to that date.
# NOTE on his: the input must contain an &his block (HIS_prefix). TokyoBay1.in includes one.
