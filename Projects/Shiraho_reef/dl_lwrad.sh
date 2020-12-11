#!/bin/bash

read -p "Enter year (YYYY): " YYYY

for ((i=1; i <= 12; i++)); do
    MM=$(printf "%02d\n" "${i}")
    wget -nc https://www.data.jma.go.jp/gmd/env/data/radiation/data/geppo/${YYYY}${MM}/DL${YYYY}${MM}_ish.txt
done
