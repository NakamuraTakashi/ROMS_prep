#!/bin/bash
# XRAIN data download and extract script
#
export TIME_START="2025-11-05T00:00"
export TIME_END="2025-11-05T23:59"
export DATA_TYPE="composite_cx"
export AREA="TOHOKU"
export OUTPUT_PATH="/mnt/f/DATA"

#
echo "Create XRAIN download data list."
uv run xrain-ls.py -u nakamura.t.av@m.titech.ac.jp -f $TIME_START -t $TIME_END $DATA_TYPE/$AREA > list.txt
echo "***** XRAIN download data list *****"
cat "./list.txt"
#
echo "Download XRAIN data."
uv run xrain-dl.py -u nakamura.t.av@m.titech.ac.jp -o $OUTPUT_PATH/xrain.tar `cat list.txt`
rm list.txt
#
cd $OUTPUT_PATH
#
echo "Unzip XRAIN data."
tar -xf xrain.tar
find . -name "*.gz" -exec gunzip -f {} +
#
rm xrain.tar
echo "XRAIN data download and extraction completed."
