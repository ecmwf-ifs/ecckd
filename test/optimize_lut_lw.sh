#!/bin/bash

set -ex

#unset OMP_NUM_THREADS

. set_paths.sh

APPLICATION=global-nwp
APP=nwp

OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_weight=0.05 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8"

TOLERANCE="0.04 0.02 0.01 0.005"
BAND_STRUCTURE="fsck wide narrow"

TOLERANCE="0.01"
BAND_STRUCTURE=wide

if [ "$APP" = nwp ]
then
    TRAINING=ckdmip_evaluation1_lw_fluxes_present.h5
    GASLIST="h2o o3 composite"
else
    TRAINING=MISSING
fi


for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}_tol${TOL}
	${BANNER} Optimizing CKD model: $MODEL_CODE

	INPUT=${WORK_LW_RAW_CKD_DIR}/lw_raw-ckd_${MODEL_CODE}.nc
	OUTPUT=${WORK_LW_CKD_DIR}/lw_ckd_${MODEL_CODE}.nc

	debug $OPTIMIZE_LUT \
	    append_path=${TRAINING_LW_FLUXES_DIR} \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    $OPTIONS \
	    "training_input=$TRAINING" \
	    gases="$GASLIST" \
	    | tee ${WORK_LW_CKD_DIR}/lw_ckd_${MODEL_CODE}.log
    done
done
