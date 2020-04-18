#!/bin/bash
# Optimize CKD look-up tables. The input requirements are the same as
# find_g_points_lw.sh

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_weight=0.05 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8"

if [ "$APP" = nwp ]
then
    TRAINING=ckdmip_evaluation1_lw_fluxes_present.h5
    GASLIST="h2o o3 composite"
else
    TRAINING=MISSING
fi

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do
    if [ "$BANDSTRUCT" = wide ]
    then
	# Map from narrow to wide bands
	BANDMAPPING="band_mapping=0 0 1 1 1 2 2 2 3 3 3 4 4"
    elif [ "$BANDSTRUCT" = fsck ]
    then
	BANDMAPPING="band_mapping=0 0 0 0 0 0 0 0 0 0 0 0 0"
    else
	BANDMAPPING="band_mapping=0 1 2 3 4 5 6 7 8 9 10 11 12"
    fi

    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}_tol${TOL}${MODEL_CODE_SUFFIX}
	${BANNER} Optimizing CKD model: $MODEL_CODE

	INPUT=${WORK_LW_RAW_CKD_DIR}/lw_raw-ckd_${MODEL_CODE}.nc
	OUTPUT=${WORK_LW_CKD_DIR}/lw_ckd_${MODEL_CODE}.nc

	$OPTIMIZE_LUT \
	    append_path=${TRAINING_LW_FLUXES_DIR} \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    $OPTIONS "$BANDMAPPING" \
	    "training_input=$TRAINING" \
	    gases="$GASLIST" \
	    | tee ${WORK_LW_CKD_DIR}/lw_ckd_${MODEL_CODE}.log
    done
done
