#!/bin/bash

DIR_NAME="2014"
LOCAL_DIR="/mnt/x/JMA_DATA/JMA_GPV"

REMOTE_DIR="/arch/jmadata/data/gpv/original"

mkdir -p ${LOCAL_DIR}/${DIR_NAME}

lftp <<EOF
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open http://database.rish.kyoto-u.ac.jp
mirror --only-missing --include-glob *_CWM_*  ${REMOTE_DIR}/${DIR_NAME} ${LOCAL_DIR}/${DIR_NAME}
EOF
