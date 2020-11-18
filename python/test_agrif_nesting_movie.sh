#!/bin/bash

################################
#SBATCH -N 1
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --constraint=HSW24
#SBATCH -J MKMOV
#SBATCH -o out_MKMOV_%J.out
#SBATCH -e err_MKMOV_%J.err
#SBATCH --time=00:55:00
####SBATCH --exclusive
################################

YEAR=2016

EXE1="${HOME}/DEV/climporn/python/nemo_movie_vect.py"
EXE2="${HOME}/DEV/climporn/python/mk_image_montage_agrif.py"

host=`hostname | cut -d '.' -f2`

if [ "${host}" = "merlat" ]; then
    CONF_I_DIR="/MEDIA/data/TROPICO05/TROPICO05-I"
    SAVE_DIR="/MEDIA/tmp/nemo"
    WRK_DIR=`pwd`
    ct0="${YEAR}1222"
    #
elif [ "${host}" = "occigen" ]; then
    CONF_I_DIR="/store/CT1/hmg2840/lbrodeau/TROPICO05/TROPICO05-I"
    SAVE_DIR="/scratch/cnt0024/hmg2840/lbrodeau/NEMO/TROPICO05/TROPICO05_NST-TRPC5N00-S/00000001-00008784"
    module purge
    module load intel/17.0 python/3.5.3
    WRK_DIR="/store/CT1/hmg2840/lbrodeau/tmp"
    ct0="${YEAR}0101"
else
    echo "PROBLEM: unknown host ==> ${host} !"; exit
fi

TTAG="6h_${YEAR}0101_${YEAR}1231"

PCONF="TROPICO05_NST" ; # parent config
CCONF="CALEDO10"      ; # child config

cdt="6"

var2plot="CURLOF"


cd ${WRK_DIR}


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


#wait
#sleep 1
#wait

echo; echo; echo; echo " TIME FOR MONTAGE !!!"; echo
${EXE2}
