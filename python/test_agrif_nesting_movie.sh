#!/bin/bash

EXE="./nemo_movie_vect.py"

CONF_I_DIR="/MEDIA/data/TROPICO05/TROPICO05-I"
SAVE_DIR="/MEDIA/tmp/nemo"
TTAG="6h_20160101_20161231"

PCONF="TROPICO05_NST" ; # parent config
CCONF="CALEDO10"      ; # child config

ct0="20161222"
cdt="6"

var2plot="CURLOF"

${EXE} -u ${SAVE_DIR}/NST/${CCONF}_${TTAG}_gridU.nc \
       -v ${SAVE_DIR}/NST/${CCONF}_${TTAG}_gridV.nc \
       -x uos -y vos -w ${var2plot} \
       -m ${CONF_I_DIR}/NST/mesh_mask_${CCONF}_trunk.nc -C ${CCONF} \
       -t ${cdt} -s ${ct0}


${EXE} -u ${SAVE_DIR}/${PCONF}-TRPC5N00_${TTAG}_gridU.nc \
       -v ${SAVE_DIR}/${PCONF}-TRPC5N00_${TTAG}_gridV.nc \
       -x uos -y vos -w ${var2plot} \
       -m ${CONF_I_DIR}/mesh_mask_TROPICO05_L31_trunk.nc -C ${PCONF} \
       -t ${cdt} -s ${ct0}

echo; echo; echo; echo " TIME FOR MONTAGE !!!"; echo
./mk_image_montage_agrif.py

