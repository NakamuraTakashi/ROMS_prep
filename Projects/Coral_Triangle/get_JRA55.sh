#!/bin/bash

read -p 'Please input your JRA55 ID: ' JRA55_ID
read -p 'Please input your Password: ' JRA55_PW

DIR_NAME="201709"
LOCAL_DIR1="Data/JRA-55/Hist/Daily/fcst_surf125"
LOCAL_DIR2="Data/JRA-55/Hist/Daily/fcst_phy2m125"

REMOTE_DIR1="JRA-55/Hist/Daily/fcst_surf125"
REMOTE_DIR2="JRA-55/Hist/Daily/fcst_phy2m125"

mkdir -p ${LOCAL_DIR1}/${DIR_NAME}
mkdir -p ${LOCAL_DIR2}/${DIR_NAME}

lftp <<END
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open -u ${JRA55_ID},${JRA55_PW} ds.data.jma.go.jp
cd ${REMOTE_DIR1}/${DIR_NAME}/
lcd ${LOCAL_DIR1}/${DIR_NAME}/
mget *
cd ../../../../../
cd ${REMOTE_DIR2}/${DIR_NAME}/
lcd ${LOCAL_DIR2}/${DIR_NAME}/
mget *
