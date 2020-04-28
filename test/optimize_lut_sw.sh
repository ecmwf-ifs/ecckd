#!/bin/bash
# Optimize CKD look-up tables. The input requirements are the same as
# find_g_points_sw.sh

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_weight=0.05 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8"

if [ "$APP" = nwp ]
then
    TRAINING=ckdmip_evaluation1_sw_fluxes_present.h5
    GASLIST="h2o o3 composite"
    INDIR=${WORK_SW_RAW_CKD_DIR}
    OUTDIR=${WORK_SW_CKD_DIR}
    INCODE=raw-ckd
    OUTCODE=ckd
elif [ "$APP" = climate ]
then
    # First pass
    TRAINING="ckdmip_evaluation1_sw_fluxes_5gas-180.h5
ckdmip_evaluation1_sw_fluxes_5gas-280.h5
ckdmip_evaluation1_sw_fluxes_5gas-415.h5
ckdmip_evaluation1_sw_fluxes_5gas-560.h5
ckdmip_evaluation1_sw_fluxes_5gas-1120.h5
ckdmip_evaluation1_sw_fluxes_5gas-2240.h5"
    GASLIST="o2n2 h2o o3 co2"
    INDIR=${WORK_SW_RAW_CKD_DIR}
    OUTDIR=${WORK_SW_RAW_CKD_DIR}
    INCODE=raw-ckd
    OUTCODE=raw2-ckd
elif [ "$APP" = climate2 ]
then
    # Second pass
    TRAINING="ckdmip_evaluation1_sw_fluxes_present.h5
ckdmip_evaluation1_sw_fluxes_ch4-350.h5
ckdmip_evaluation1_sw_fluxes_ch4-700.h5
ckdmip_evaluation1_sw_fluxes_ch4-1200.h5
ckdmip_evaluation1_sw_fluxes_ch4-2600.h5
ckdmip_evaluation1_sw_fluxes_ch4-3500.h5
ckdmip_evaluation1_sw_fluxes_n2o-190.h5
ckdmip_evaluation1_sw_fluxes_n2o-270.h5
ckdmip_evaluation1_sw_fluxes_n2o-405.h5
ckdmip_evaluation1_sw_fluxes_n2o-540.h5
ckdmip_evaluation1_sw_fluxes_cfc11-0.h5
ckdmip_evaluation1_sw_fluxes_cfc11-2000.h5"
# Don't bother optimizing CFC12 - more accurate out of the box
#ckdmip_evaluation1_sw_fluxes_cfc12-0.h5
#ckdmip_evaluation1_sw_fluxes_cfc12-550.h5
#    EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_5gas-415.h5"
    EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_rel-415.h5"
    GASLIST="ch4 n2o cfc11"
    INDIR=${WORK_SW_RAW_CKD_DIR}
    OUTDIR=${WORK_SW_CKD_DIR}
    INCODE=raw2-ckd
    OUTCODE=ckd
elif [ "$APP" = climate3 ]
then
    # All-in-one
    GASLIST="o2n2 h2o o3 co2 ch4 n2o cfc11"
    INDIR=${WORK_SW_RAW_CKD_DIR}
    OUTDIR=${WORK_SW_CKD_DIR}
    INCODE=raw-ckd
    OUTCODE=ckd
    TRAINING="ckdmip_evaluation1_sw_fluxes_present.h5
ckdmip_evaluation1_sw_fluxes_co2-180.h5
ckdmip_evaluation1_sw_fluxes_co2-280.h5
ckdmip_evaluation1_sw_fluxes_co2-560.h5
ckdmip_evaluation1_sw_fluxes_co2-1120.h5
ckdmip_evaluation1_sw_fluxes_co2-2240.h5
ckdmip_evaluation1_sw_fluxes_ch4-350.h5
ckdmip_evaluation1_sw_fluxes_ch4-700.h5
ckdmip_evaluation1_sw_fluxes_ch4-1200.h5
ckdmip_evaluation1_sw_fluxes_ch4-2600.h5
ckdmip_evaluation1_sw_fluxes_ch4-3500.h5
ckdmip_evaluation1_sw_fluxes_n2o-190.h5
ckdmip_evaluation1_sw_fluxes_n2o-270.h5
ckdmip_evaluation1_sw_fluxes_n2o-405.h5
ckdmip_evaluation1_sw_fluxes_n2o-540.h5
ckdmip_evaluation1_sw_fluxes_cfc11-0.h5
ckdmip_evaluation1_sw_fluxes_cfc11-2000.h5"
elif [ "$APP" = climate4 ]
then
    # First pass
    TRAINING="ckdmip_evaluation1_sw_fluxes_rel-180.h5
ckdmip_evaluation1_sw_fluxes_rel-280.h5
ckdmip_evaluation1_sw_fluxes_rel-415.h5
ckdmip_evaluation1_sw_fluxes_rel-560.h5
ckdmip_evaluation1_sw_fluxes_rel-1120.h5
ckdmip_evaluation1_sw_fluxes_rel-2240.h5"
    GASLIST="composite h2o o3 co2"
#    MODEL_CODE_SUFFIX=-rel
    INDIR=${WORK_SW_RAW_CKD_DIR}
    OUTDIR=${WORK_SW_RAW_CKD_DIR}
    INCODE=raw-ckd
    OUTCODE=raw2-ckd
else 
    echo "Error"
    exit 1
fi

mkdir -p ${OUTDIR}

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

	INPUT=${INDIR}/sw_${INCODE}_${MODEL_CODE}.nc
	OUTPUT=${OUTDIR}/sw_${OUTCODE}_${MODEL_CODE}.nc
	LOG=${OUTDIR}/sw_${OUTCODE}_${MODEL_CODE}.log

	$OPTIMIZE_LUT \
	    append_path=${TRAINING_SW_FLUXES_DIR}:${WORK_SW_LBL_FLUX_DIR} \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    $OPTIONS "$BANDMAPPING" \
	    training_input="$TRAINING" \
	    gases="$GASLIST" $EXTRA_ARGS \
	    | tee $LOG
    done
done
