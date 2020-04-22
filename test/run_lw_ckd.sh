#!/bin/bash
# Run two-stream radiative transfer. The input requirements are the same as
# find_g_points_lw.sh 

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

#VERSIONS="raw-ckd ckd"
VERSIONS=ckd

if [ "$APP" = nwp ]
then
    SCENARIOS=present
else
#    VERSIONS="raw-ckd raw2-ckd ckd"
    SCENARIOS="glacialmax
preindustrial
present
future
co2-180
co2-280
co2-560
co2-1120
co2-2240
ch4-350
ch4-700
ch4-1200
ch4-2600
ch4-3500
n2o-190
n2o-270
n2o-405
n2o-540
cfc11-0
cfc11-2000
cfc12-0
cfc12-550
co2-180-ch4-350
co2-2240-ch4-350
co2-180-ch4-3500
co2-2240-ch4-3500
co2-180-n2o-190
co2-2240-n2o-190
co2-180-n2o-540
co2-2240-n2o-540
ch4-350-n2o-190
ch4-3500-n2o-190
ch4-350-n2o-540
ch4-3500-n2o-540"

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
		#OUTPUT=${WORK_LW_CKD_OD_DIR}/lw_${VER}_${MODEL_CODE}_optical-depth_${SCENARIO}.nc
		OUTPUT=${WORK_LW_CKD_OD_DIR}/ecckd-${VER}_evaluation1_lw_${APPLICATION}_${BANDSTRUCT}-tol${TOL}_optical-depth_${SCENARIO}.nc
		$LW_CKD ckd_model=$CKD_MODEL input=${INPUT} output=${OUTPUT}
	    done
	done
    done
done
