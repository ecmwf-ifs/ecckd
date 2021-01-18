#!/bin/bash
# Reorder the spectra of H2O, O3 and the composite of well-mixed gases

. config.h

mkdir -p ${WORK_SW_ORDER_DIR}

#OPTIONS="ssi=$MMM_SW_SSI threshold_optical_depth=0.5"
OPTIONS="ssi=$MMM_SW_SSI threshold_optical_depth=0.25"

GAS_LIST="composite h2o_median o3_median co2_present ch4_present n2o_present o2n2_constant rayleigh_present"

# Loop through the median/present concentrations of each gas and
# reorder
for GAS_SCENARIO in ${GAS_LIST}
do
    GAS=$(echo ${GAS_SCENARIO} | awk -F_ '{print $1}')
    GAS1=$(echo ${GAS} | awk -F- '{print $1}')

    if [ "$GAS1" = composite ]
    then
	INPUT=${WELL_MIXED_SW_SPECTRA}
    elif [ "$GAS1" = o2n2 ]
    then
	INPUT=${WELL_MIXED_SW_SPECTRA_O2N2}
    else
	INPUT=${MMM_SW_SPECTRA_DIR}/ckdmip_mmm_sw_spectra_${GAS_SCENARIO}.h5
    fi

    for BANDSTRUCT in ${BAND_STRUCTURE}
    do

	OUTPUT=${WORK_SW_ORDER_DIR}/sw_order_${BANDSTRUCT}_${GAS}.h5

	if [ ! -f ${OUTPUT} ]
	then
	    ${BANNER} Reordering ${GAS}, band structure ${BANDSTRUCT}
	    if [ "$BANDSTRUCT" = narrow ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ${OPTIONS} \
		    "wavenumber1=$WN1_SW_NARROW" "wavenumber2=$WN2_SW_NARROW"
	    elif [ "$BANDSTRUCT" = wide ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ${OPTIONS} \
		    "wavenumber1=$WN1_SW_WIDE" "wavenumber2=$WN2_SW_WIDE"
	    elif [ "$BANDSTRUCT" = double ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ${OPTIONS} \
		    "wavenumber1=$WN1_SW_DOUBLE" "wavenumber2=$WN2_SW_DOUBLE"
	    elif [ "$BANDSTRUCT" = rgb ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ${OPTIONS} \
		    "wavenumber1=$WN1_SW_RGB" "wavenumber2=$WN2_SW_RGB"
	    elif [ "$BANDSTRUCT" = fsck ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ${OPTIONS}
	    else
		${BANNER_ERROR} "Band structure\"$BANDSTRUCT\" not understood"
		exit 1
	    fi
	else
	    ${BANNER_SKIP} Skipping reordering of ${GAS} as file present: ${OUTPUT}
	fi

    done
done
