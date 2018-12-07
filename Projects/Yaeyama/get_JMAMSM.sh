#!/bin/bash

DIR_NAME="2017"
MONTH="05"

LOCAL_DIR1="N:/Data/JMA-MSM/MSM-S"
LOCAL_DIR2="N:/Data/JMA-MSM/MSM-S/r1h"

REMOTE_DIR1="MSM-S"
REMOTE_DIR2="MSM-S/r1h"

mkdir ${LOCAL_DIR1}/${DIR_NAME}
mkdir ${LOCAL_DIR2}/${DIR_NAME}

lftp <<END
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/netcdf/
cd ${REMOTE_DIR1}/${DIR_NAME}/
lcd ${LOCAL_DIR1}/${DIR_NAME}/
mget ${MONTH}*
cd ../../
cd ${REMOTE_DIR2}/${DIR_NAME}/
lcd ${LOCAL_DIR2}/${DIR_NAME}/
mget ${MONTH}*
