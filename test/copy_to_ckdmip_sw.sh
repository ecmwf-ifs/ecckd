#!/bin/bash

. config.h

# Get ncdump if on ECMWF system
if [ "$ECACCOUNT" ]
then
    module load netcdf4
fi

VERSIONS="ckd raw"
APPLICATION=climate
APPLICATION=global-nwp
#APPLICATION=limited-area-nwp
BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.16 0.08 0.04 0.02 0.01 0.005"

BAND_STRUCTURE=fsck
TOLERANCE=4

BAND_STRUCTURE=wide
TOLERANCE=0.1

mkdir -p ${CKDMIP_RESULTS_DIR}/sw_spectral-definition/
mkdir -p ${CKDMIP_RESULTS_DIR}/sw_optical-depth/
mkdir -p ${CKDMIP_RESULTS_DIR}/sw_fluxes/

# Loop over each band structure, tolerance, version (raw or final) and
# scenario
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	for VER in $VERSIONS
	do

	    if [ ! "$VER" = ckd ]
	    then
		VER_SUFFIX=-$VER
	    else
		unset VER_SUFFIX
	    fi

	    MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}

	    echo From directory ${WORK_SW_CKD_DIR}
	    cd ${WORK_SW_CKD_DIR}

	    CKD_FILE=${ECCKD_PREFIX}_sw_ckd-definition_${MODEL_CODE}.nc

	    # Get number of g points
	    NG=$(ncdump -h $CKD_FILE | head -10 | grep g_point | awk '{print $3}')

	    NEW_MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-${NG}

	    NEW_CKD_FILE=${CKDMIP_RESULTS_DIR}/sw_spectral-definition/${ECCKD_PREFIX}_sw_${NEW_MODEL_CODE}${VER_SUFFIX}_spectral-definition.nc

	    echo "  Copying $CKD_FILE -> $NEW_CKD_FILE"
	    cp -f $CKD_FILE $NEW_CKD_FILE

	    for FILE_TYPE in optical-depth fluxes
	    do
		if [ "$FILE_TYPE" = optical-depth ]
		then
		    echo "From directory ${WORK_SW_CKD_OD_DIR}"
		    cd ${WORK_SW_CKD_OD_DIR}
		else
		    echo "From directory ${WORK_SW_FLUX_DIR}"
		    cd ${WORK_SW_FLUX_DIR}
		fi

		for FILE in ${ECCKD_PREFIX}_${EVALUATION_CODE}_sw_${MODEL_CODE}${VER_SUFFIX}_${FILE_TYPE}_*.nc
		do
		    NEW_FILE=$(echo $FILE | sed "s|${MODEL_CODE}|${NEW_MODEL_CODE}|")
		    echo "  Copying $FILE -> ${CKDMIP_RESULTS_DIR}/sw_${FILE_TYPE}/$NEW_FILE"
		    cp -f $FILE ${CKDMIP_RESULTS_DIR}/sw_${FILE_TYPE}/$NEW_FILE
		done

	    done

	done
    done
done
