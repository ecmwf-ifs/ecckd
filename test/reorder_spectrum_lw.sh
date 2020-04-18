#!/bin/bash
# Reorder the spectra of H2O, O3 and the composite of well-mixed gases

. config.h

mkdir -p ${WORK_LW_ORDER_DIR}

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
		    "wavenumber1=$WN1_LW_NARROW" "wavenumber2=$WN2_LW_NARROW"
	    elif [ "$BAND_STRUCTURE" = wide ]
	    then
		${REORDER_SPECTRUM_LW} iprofile=0 input=$INPUT output=$OUTPUT \
		    "wavenumber1=$WN1_LW_WIDE" "wavenumber2=$WN2_LW_WIDE"
	    else
		# Assuming FSCK
		${REORDER_SPECTRUM_LW} iprofile=0 input=$INPUT output=$OUTPUT
	    fi
	else
	    ${BANNER_SKIP} Skipping reordering of ${GAS} as file present: ${OUTPUT}
	fi

    done
done
