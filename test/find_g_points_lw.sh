#!/bin/bash

unset OMP_NUM_THREADS
#export OMP_NUM_THREADS=6
set -e

. set_paths.sh

# Interpret first argument as a script to be sourced
if [ "$#" -gt 0 ]
then
    . $1
fi

if [ ! ${MIN_PRESSURE} ]
then
    ${BANNER_ERROR} '"MIN_PRESSURE"' not specified
    exit 1
fi

if [ ! "${TOLERANCE}" ]
then
    ${BANNER_ERROR} '"TOLERANCE"' not specified
    exit 1
fi

if [ ! ${APPLICATION} ]
then
    ${BANNER_ERROR} '"APPLICATION"' not specified
    exit 1
fi

if [ ! ${APP} ]
then
    ${BANNER_ERROR} '"APP"' not specified
    exit 1
fi

if [ ! "${BAND_STRUCTURE}" ]
then
    ${BANNER_ERROR} '"BAND_STRUCTURE"' not specified
    exit 1
fi

mkdir -p ${WORK_LW_GPOINTS_DIR}

for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do

	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}_tol${TOL}

	${BANNER} Finding g-points: $MODEL_CODE

debug	${FIND_G_POINTS_LW} \
	    append_path="${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}" \
	    heating_rate_tolerance=${TOL} \
	    min_pressure=${MIN_PRESSURE} \
	    h2o.reordering_input=lw_order_${BANDSTRUCT}_h2o.h5 \
	    o3.reordering_input=lw_order_${BANDSTRUCT}_o3.h5 \
	    composite.reordering_input=lw_order_${BANDSTRUCT}_composite.h5 \
	    output=${WORK_LW_GPOINTS_DIR}/lw_gpoints_${MODEL_CODE}.h5 \
	    config_find_g_points_lw_${APP}.cfg \
	    | tee ${WORK_LW_GPOINTS_DIR}/lw_gpoints_${MODEL_CODE}.log
    done
done
