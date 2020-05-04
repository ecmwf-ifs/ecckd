#!/bin/bash

. config.h

# Get ncdump if on ECMWF system
if [ "$ECACCOUNT" ]
then
    module load netcdf4
fi

VERSIONS="ckd"
APPLICATION=climate
APPLICATION=global-nwp
APPLICATION=limited-area-nwp
BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.16 0.08 0.04 0.02 0.01 0.005"

mkdir -p ${CKDMIP_RESULTS_DIR}/lw_spectral-definition/
mkdir -p ${CKDMIP_RESULTS_DIR}/lw_optical-depth/

# Loop over each band structure, tolerance, version (raw or final) and
# scenario
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	for VER in $VERSIONS
	do
	    MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}

	    echo From directory ${WORK_LW_CKD_DIR}
	    cd ${WORK_LW_CKD_DIR}

	    CKD_FILE=${ECCKD_PREFIX}_lw_ckd-definition_${MODEL_CODE}.nc

	    # Get number of g points
	    NG=$(ncdump -h $CKD_FILE | head -10 | grep g_point | awk '{print $3}')

	    NEW_MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-${NG}

	    NEW_CKD_FILE=${CKDMIP_RESULTS_DIR}/lw_spectral-definition/${ECCKD_PREFIX}_lw_${NEW_MODEL_CODE}_spectral-definition.nc

	    echo "  Copying $CKD_FILE -> $NEW_CKD_FILE"
	    cp -f $CKD_FILE $NEW_CKD_FILE

	    echo "From directory ${WORK_LW_CKD_OD_DIR}"
	    cd ${WORK_LW_CKD_OD_DIR}

	    for FILE in ${ECCKD_PREFIX}_${EVALUATION_CODE}_lw_${MODEL_CODE}_optical-depth_*.nc
	    do
		NEW_FILE=$(echo $FILE | sed "s|${MODEL_CODE}|${NEW_MODEL_CODE}|")
		echo "  Copying $FILE -> ${CKDMIP_RESULTS_DIR}/lw_optical-depth/$NEW_FILE"
		cp -f $FILE ${CKDMIP_RESULTS_DIR}/lw_optical-depth/$NEW_FILE
	    done

	done
    done
done
