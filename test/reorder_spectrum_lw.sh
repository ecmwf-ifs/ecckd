#!/bin/bash

set -e

. set_paths.sh

mkdir -p ${WORK_LW_ORDER_DIR}

# Loop through the median/present concentrations of each gas and
# reorder
for GAS_SCENARIO in merge-well-mixed h2o_median o3_median
do
    GAS=$(echo ${GAS_SCENARIO} | awk -F_ '{print $1}')
    GAS1=$(echo ${GAS} | awk -F- '{print $1}')

    if [ "$GAS1" = merge ]
    then
	INPUT=${WELL_MIXED_LW_SPECTRA}
    else
	INPUT=${MMM_LW_SPECTRA_DIR}/ckdmip_mmm_lw_spectra_${GAS_SCENARIO}.h5
    fi
    OUTPUT=${WORK_LW_ORDER_DIR}/lw_order_${GAS}.h5

    if [ ! -f ${OUTPUT} ]
    then
	${BANNER} Reordering ${GAS}
	${REORDER_SPECTRUM_LW} iprofile=0 input=$INPUT output=$OUTPUT
    else
	${BANNER_SKIP} Skipping reordering of ${GAS} as file present: ${OUTPUT}
    fi
done
