#!/bin/bash

DIR_NAME="2025"
#LOCAL_DIR="//gshs.t4.gsic.titech.ac.jp/t4_bs/tga-NakamuLab/DATA/JMA_DATA/JMA_GPV"
LOCAL_DIR="/mnt/x/JMA_DATA/JMA_GPV"

REMOTE_DIR="/arch/jmadata/data/gpv/original"

mkdir -p ${LOCAL_DIR}/${DIR_NAME}

lftp <<EOF
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open http://database.rish.kyoto-u.ac.jp
mirror --only-missing --include-glob *_CWM_* --include-glob *_MSM_GPV_Rjp_Lsurf_FH00-15_* --include-glob MSM*SFC018_* ${REMOTE_DIR}/${DIR_NAME} ${LOCAL_DIR}/${DIR_NAME}
EOF
