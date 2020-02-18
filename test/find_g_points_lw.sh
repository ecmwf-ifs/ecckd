#!/bin/bash

set -e

. set_paths.sh

${BANNER} Finding g-points

if [ ! ${MODEL_CODE} ]
then
    ${BANNER_ERROR} '"MODEL_CODE"' not specified
    exit 1
fi

if [ ! ${MIN_PRESSURE} ]
then
    ${BANNER_ERROR} '"MIN_PRESSURE"' not specified
    exit 1
fi

if [ ! ${TOLERANCE} ]
then
    ${BANNER_ERROR} '"TOLERANCE"' not specified
    exit 1
fi

if [ ! ${APP} ]
then
    ${BANNER_ERROR} '"APP"' not specified
    exit 1
fi

mkdir -p ${WORK_LW_GPOINTS_DIR}

debug ${FIND_G_POINTS_LW} \
    append_path="${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}" \
    heating_rate_tolerance=${TOLERANCE} \
    min_pressure=${MIN_PRESSURE} repartition_factor=1 repartition_repeat=5 tolerance_tolerance=0.05 \
    output=${WORK_LW_GPOINTS_DIR}/lw_gpoints_${MODEL_CODE}.h5 \
    config_find_g_points_lw_${APP}.cfg
