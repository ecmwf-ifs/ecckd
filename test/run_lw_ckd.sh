#!/bin/bash
# Run two-stream radiative transfer. The input requirements are the same as
# find_g_points_lw.sh 

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

VERSIONS="raw-ckd ckd"
#VERSIONS="ckd"
#OUTPUT_SUFFIX=-t2

if [ "$APP" = nwp ]
then
    SCENARIOS=present
else
    SCENARIOS=
fi

# Loop over each band structure, tolerance, version (raw or final) and
# scenario
for BANDSTRUCT in $BAND_STRUCTURE
do

    for TOL in $TOLERANCE
    do
	for VER in $VERSIONS
	do
	    for SCENARIO in $SCENARIOS
	    do
		MODEL_CODE=${APPLICATION}_${BANDSTRUCT}_tol${TOL}${MODEL_CODE_SUFFIX}
		CKD_MODEL=${WORK_DIR}/lw_${VER}/lw_${VER}_${MODEL_CODE}.nc
		INPUT=${TRAINING_CONC_DIR}/ckdmip_${TRAINING_CODE}_concentrations_${SCENARIO}.nc
		OUTPUT=${WORK_LW_CKD_OD_DIR}/lw_${VER}_${MODEL_CODE}_optical-depth_${SCENARIO}.nc
		$LW_CKD ckd_model=$CKD_MODEL input=${INPUT} output=${OUTPUT}
	    done
	done
    done
done
