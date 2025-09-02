#!/bin/bash
# Please set your username and password
# **CAUTION: Do not upload your username and password publicly.**
#    Please rename this script file to "my_TE-download.sh" 
#    after setting your username and password.
#    Then change the permission of this file to -rwx------ (700).
#    $ chmod 700 my_TE-download.sh
USERNAME="******"
PASSWORD="******"

DIR_NAME="2025"
LOCAL_DIR="/mnt/f/DATA/TodaysEarth/TE-japan/MSM/hourly"
REMOTE_DIR="/TE-japan/MSM/hourly"

mkdir -p ${LOCAL_DIR}/${DIR_NAME}

lftp <<EOF
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open -u ${USERNAME},${PASSWORD} ftp.eorc.jaxa.jp
mirror --only-missing --parallel=4 \
--include-glob *RIVOUT.nc \
--include-glob *GPRCT.nc \
--include-glob *GSNWL.nc \
--include-glob *SSRD.nc \
--include-glob *SLRD.nc \
${REMOTE_DIR}/${DIR_NAME} ${LOCAL_DIR}/${DIR_NAME}
EOF

