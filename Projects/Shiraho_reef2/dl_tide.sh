#!/bin/bash

read -p "Enter year (YYYY): " YYYY

for ((i=1; i <= 12; i++)); do
    MM=$(printf "%02d\n" "${i}")
    if [ ${YYYY} -gt 1997 ]; then
        wget -nc https://www.data.jma.go.jp/gmd/kaiyou/data/db/tide/genbo/${YYYY}/${YYYY}${MM}/hry${YYYY}${MM}IS.txt
    else
        if [ ${YYYY} = 1997 ] && [ ${i} -ge 4 ]; then
            wget -nc https://www.data.jma.go.jp/gmd/kaiyou/data/db/tide/genbo/${YYYY}/${YYYY}${MM}/hry${YYYY}${MM}IS.txt
        else
            wget -nc https://www.data.jma.go.jp/gmd/kaiyou/data/db/tide/sea_lev_var/${YYYY}/${YYYY}${MM}/hry${YYYY}${MM}IS.txt
        fi
    fi
done
