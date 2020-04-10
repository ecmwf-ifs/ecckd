#!/bin/bash

set -e

unset OMP_NUM_THREADS

. set_paths.sh

mkdir -p ${WORK_LW_ORDER_DIR}

WN1_NARROW="0 350 500 630 700 820 980 1080 1180 1390 1480 1800 2080"
WN2_NARROW="350 500 630 700 820 980 1080 1180 1390 1480 1800 2080 3260"
WN1_WIDE="0 500 820 1180 1800"
WN2_WIDE="500 820 1180 1800 3260"


# Loop through the median/present concentrations of each gas and
# reorder
for GAS_SCENARIO in composite h2o_median o3_median
do
    GAS=$(echo ${GAS_SCENARIO} | awk -F_ '{print $1}')
    GAS1=$(echo ${GAS} | awk -F- '{print $1}')

    if [ "$GAS1" = composite ]
    then
	INPUT=${WELL_MIXED_LW_SPECTRA}
    else
	INPUT=${MMM_LW_SPECTRA_DIR}/ckdmip_mmm_lw_spectra_${GAS_SCENARIO}.h5
    fi

    for BAND_STRUCTURE in fsck wide narrow
    do

	OUTPUT=${WORK_LW_ORDER_DIR}/lw_order_${BAND_STRUCTURE}_${GAS}.h5

	if [ ! -f ${OUTPUT} ]
	then
	    ${BANNER} Reordering ${GAS}, band structure ${BAND_STRUCTURE}
	    if [ "$BAND_STRUCTURE" = narrow ]
	    then
		${REORDER_SPECTRUM_LW} iprofile=0 input=$INPUT output=$OUTPUT \
		    "wavenumber1=$WN1_NARROW" "wavenumber2=$WN2_NARROW"
	    elif [ "$BAND_STRUCTURE" = wide ]
	    then
		${REORDER_SPECTRUM_LW} iprofile=0 input=$INPUT output=$OUTPUT \
		    "wavenumber1=$WN1_WIDE" "wavenumber2=$WN2_WIDE"
	    else
		# Assuming FSCK
		${REORDER_SPECTRUM_LW} iprofile=0 input=$INPUT output=$OUTPUT
	    fi
	else
	    ${BANNER_SKIP} Skipping reordering of ${GAS} as file present: ${OUTPUT}
	fi

    done
done
