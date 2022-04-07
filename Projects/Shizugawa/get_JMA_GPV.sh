#!/bin/bash

DIR_NAME="2021"
LOCAL_DIR="/cygdrive/o/JMA_DATA/JMA_GPV"

REMOTE_DIR="/arch/jmadata/data/gpv/original"

mkdir -p ${LOCAL_DIR}/${DIR_NAME}

lftp <<EOF
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open http://database.rish.kyoto-u.ac.jp
mirror --only-missing --include-glob *_CWM_* --include-glob *_MSM_GPV_Rjp_Lsurf_FH00-15_* --include-glob MSM*SFC018_* ${REMOTE_DIR}/${DIR_NAME} ${LOCAL_DIR}/${DIR_NAME}
EOF
