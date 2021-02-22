#!/bin/bash
read -p "Enter year (YYYY): " YYYY

LOCAL_DIR="/cygdrive/h/WW3"

mkdir -p ${LOCAL_DIR}

cd ${LOCAL_DIR}

for ((i=1; i <= 12; i++)); do
    MM=$(printf "%02d\n" "${i}")
    wget -nc https://polar.ncep.noaa.gov/waves/hindcasts/nopp-phase2/${YYYY}${MM}/gribs/multi_reanal.glo_30m_ext.hs.${YYYY}${MM}.grb2
    wget -nc https://polar.ncep.noaa.gov/waves/hindcasts/nopp-phase2/${YYYY}${MM}/gribs/multi_reanal.glo_30m_ext.dp.${YYYY}${MM}.grb2
    wget -nc https://polar.ncep.noaa.gov/waves/hindcasts/nopp-phase2/${YYYY}${MM}/gribs/multi_reanal.glo_30m_ext.tp.${YYYY}${MM}.grb2
done

