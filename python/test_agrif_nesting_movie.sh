#!/bin/bash

EXE1="${HOME}/DEV/climporn/python/nemo_movie_vect.py"
EXE2="${HOME}/DEV/climporn/python/mk_image_montage_agrif.py"


host=`hostname | cut -d '.' -f2`

if [ "${host}" = "merlat" ]; then
    CONF_I_DIR="/MEDIA/data/TROPICO05/TROPICO05-I"
    SAVE_DIR="/MEDIA/tmp/nemo"
elif [ "${host}" = "occigen" ]; then
    CONF_I_DIR="/store/CT1/hmg2840/lbrodeau/TROPICO05/TROPICO05-I"
    SAVE_DIR="/scratch/cnt0024/hmg2840/lbrodeau/NEMO/TROPICO05/TROPICO05_NST-TRPC5N00-S/00000001-00008784"
else
    echo "PROBLEM: unknown host ==> ${host} !"; exit
fi

    
TTAG="6h_20160101_20161231"

PCONF="TROPICO05_NST" ; # parent config
CCONF="CALEDO10"      ; # child config

ct0="20161222"
cdt="6"

var2plot="CURLOF"

${EXE1} -u ${SAVE_DIR}/NST/${CCONF}_${TTAG}_gridU.nc \
       -v ${SAVE_DIR}/NST/${CCONF}_${TTAG}_gridV.nc \
       -x uos -y vos -w ${var2plot} \
       -m ${CONF_I_DIR}/NST/mesh_mask_${CCONF}_trunk.nc -C ${CCONF} \
       -t ${cdt} -s ${ct0}


${EXE1} -u ${SAVE_DIR}/${PCONF}-TRPC5N00_${TTAG}_gridU.nc \
       -v ${SAVE_DIR}/${PCONF}-TRPC5N00_${TTAG}_gridV.nc \
       -x uos -y vos -w ${var2plot} \
       -m ${CONF_I_DIR}/mesh_mask_TROPICO05_L31_trunk.nc -C ${PCONF} \
       -t ${cdt} -s ${ct0}

echo; echo; echo; echo " TIME FOR MONTAGE !!!"; echo
${EXE2}

