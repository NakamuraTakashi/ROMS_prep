#!/bin/bash
#
# TPXO9-ATLAS
# http://volkov.oce.orst.edu/tides/tpxo9_atlas.html

lftp <<END
set ftp:proxy http://proxy.noc.titech.ac.jp:3128/
open ftp://ftp.oce.orst.edu/dist/tides/TPXO9_atlas/
mget *_tpxo9_atlas_30
