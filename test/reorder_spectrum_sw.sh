#!/bin/bash
# Reorder the spectra of H2O, O3 and the composite of well-mixed gases

. config.h

mkdir -p ${WORK_SW_ORDER_DIR}

# Loop through the median/present concentrations of each gas and
# reorder
for GAS_SCENARIO in composite h2o_median o3_median
do
    GAS=$(echo ${GAS_SCENARIO} | awk -F_ '{print $1}')
    GAS1=$(echo ${GAS} | awk -F- '{print $1}')

    if [ "$GAS1" = composite ]
    then
	INPUT=${WELL_MIXED_SW_SPECTRA}
    else
	INPUT=${MMM_SW_SPECTRA_DIR}/ckdmip_mmm_sw_spectra_${GAS_SCENARIO}.h5
    fi

    for BAND_STRUCTURE in fsck wide narrow
    do

	OUTPUT=${WORK_SW_ORDER_DIR}/sw_order_${BAND_STRUCTURE}_${GAS}.h5

	if [ ! -f ${OUTPUT} ]
	then
	    ${BANNER} Reordering ${GAS}, band structure ${BAND_STRUCTURE}
	    if [ "$BAND_STRUCTURE" = narrow ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ssi=$MMM_SW_SSI \
		    "wavenumber1=$WN1_SW_NARROW" "wavenumber2=$WN2_SW_NARROW"
	    elif [ "$BAND_STRUCTURE" = wide ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ssi=$MMM_SW_SSI \
		    "wavenumber1=$WN1_SW_WIDE" "wavenumber2=$WN2_SW_WIDE"
	    else
		# Assuming FSCK
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ssi=$MMM_SW_SSI
	    fi
	else
	    ${BANNER_SKIP} Skipping reordering of ${GAS} as file present: ${OUTPUT}
	fi

    done
done