#!/bin/bash
#
# (C) Copyright 2019- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

. config.h

# Get ncdump if on ECMWF system
if [ "$ECACCOUNT" ]
then
    module load netcdf4
fi

mkdir -p ${CKDMIP_RESULTS_DIR}/sw_spectral-definition/
mkdir -p ${CKDMIP_RESULTS_DIR}/sw_optical-depth/
mkdir -p ${CKDMIP_RESULTS_DIR}/sw_fluxes/

# If the version of the model to run is not stated, do just the final
# one
if [ -z "$VERSIONS" ]
then
    VERSIONS=ckd
fi

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

	    if [ ! -r "$CKD_FILE" ]
	    then
		echo ... not found, trying ${WORK_SW_RAW_CKD_DIR}
		cd ${WORK_SW_RAW_CKD_DIR}
		CKD_FILE=${ECCKD_PREFIX}_sw_${VER}-ckd-definition_${MODEL_CODE}.nc
	    fi

	    # Get number of g points
	    NG=$(ncdump -h $CKD_FILE | head -10 | grep g_point | awk '{print $3}')

	    # Use "b" suffix to indicate that training used both
	    # evaluation1 and evaluation2 LBL reference data
	    if [ "$TRAINING_BOTH" = yes ]
	    then
		NG=${NG}b
	    fi

	    NEW_MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-${NG}${MODEL_CODE_SUFFIX}

	    NEW_CKD_FILE=${CKDMIP_RESULTS_DIR}/sw_spectral-definition/${ECCKD_PREFIX}_sw_${NEW_MODEL_CODE}_spectral-definition.nc

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
