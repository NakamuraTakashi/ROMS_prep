#!/bin/bash
# ============================================================================================
#  Interpolate ROMS atmospheric forcing from the TokyoBay2 grid onto the nested
#  TokyoBay3 grid (horizontal bilinear remapping with CDO).
#
#  Input  : TokyoBay2_frc_MSMgb_*.nc   (on TokyoBay2 rho grid, 192 x 137)
#  Output : TokyoBay3_frc_MSMgb_*.nc   (on TokyoBay3 rho grid, 269 x 293)
#
#  Notes
#   - ROMS forcing files store only (time, eta_rho, xi_rho); the lon/lat live in the GRID
#     files, so we attach them to make CDO see a curvilinear grid (cdo setgrid).
#   - Atmospheric fields are valid everywhere (no land mask), and both grids have angle = 0
#     (aligned to east/north), so Uwind/Vwind need no rotation -> plain bilinear is correct.
#   - TokyoBay3 is fully nested inside TokyoBay2, so every target point is filled (no gaps).
#
#  Requires: cdo, ncks, ncatted (NCO).
# ============================================================================================
set -eu

# ----- Paths (edit if needed) -----
DATADIR=../../../COAWST_DATA/Shizugawa
SRC_GRID=${DATADIR}/Shizugawa1/Grid/Shizugawa1_grd_v0.1.nc
DST_GRID=${DATADIR}/Shizugawa3/Grid/Shizugawa3_grd_v0.5.nc

INPUT_GLOB="Shizugawa1_frc_MSMgb_*.nc"   # forcing files to convert (on TokyoBay2 grid)
SRC_TAG="Shizugawa1"                       # replaced by DST_TAG to build output names
DST_TAG="Shizugawa3"
OUTDIR="."                                # where to write the TokyoBay3 forcing files

REMAP=remapbil                            # remapbil (bilinear); use remapnn for nearest-neighbour

# ----- Checks -----
for cmd in cdo ncks ncatted; do command -v "$cmd" >/dev/null || { echo "ERROR: '$cmd' not found"; exit 1; }; done
for f in "$SRC_GRID" "$DST_GRID"; do [ -f "$f" ] || { echo "ERROR: grid not found: $f"; exit 1; }; done

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ----- 1) Build CF grid-description files so CDO reads the ROMS rho grids as curvilinear -----
#         (carry lon_rho/lat_rho + one variable that points to them via 'coordinates')
build_grid() {  # $1 = ROMS grid file, $2 = output CF grid file
  ncks -O -v lon_rho,lat_rho,h "$1" "$2"
  ncatted -O -a coordinates,h,c,c,"lat_rho lon_rho" "$2"
}
build_grid "$SRC_GRID" "$TMP/src_grid.nc"
build_grid "$DST_GRID" "$TMP/dst_grid.nc"

# ----- 2) Generate remap weights once (same src->dst grids for every file) -----
shopt -s nullglob
files=( $INPUT_GLOB )
[ ${#files[@]} -gt 0 ] || { echo "ERROR: no input files match '$INPUT_GLOB'"; exit 1; }

gen_op=${REMAP/remap/gen}                  # remapbil->genbil, remapnn->gennn
echo "Generating ${REMAP} weights (${gen_op}) ..."
cdo -s "${gen_op},${TMP}/dst_grid.nc" -setgrid,"${TMP}/src_grid.nc" -seltimestep,1 "${files[0]}" "$TMP/weights.nc"

# ----- 3) Remap each forcing file, then restore the plain ROMS-forcing layout -----
mkdir -p "$OUTDIR"
for in in "${files[@]}"; do
  out="${OUTDIR}/$(basename "${in/$SRC_TAG/$DST_TAG}")"
  echo "Remapping ${in}  ->  ${out}"
  cdo -s remap,"${TMP}/dst_grid.nc","${TMP}/weights.nc" -setgrid,"${TMP}/src_grid.nc" "$in" "$TMP/remap.nc"
  # CDO adds lon_rho/lat_rho coordinate variables + a 'coordinates' attribute; drop them so
  # the output matches the original ROMS forcing files (time, eta_rho, xi_rho only).
  # (To KEEP the lon/lat in the output, replace the next two lines with: cp "$TMP/remap.nc" "$out")
  ncks -O -C -x -v lon_rho,lat_rho "$TMP/remap.nc" "$out"
  ncatted -O -a coordinates,,d,, "$out"
done

echo "Done. Output written to: ${OUTDIR}"
