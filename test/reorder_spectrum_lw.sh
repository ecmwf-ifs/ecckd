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
#
# This script reorders the longwave spectra of H2O, O3 and the
# composite of well-mixed gases. It is called from do_all_lw.sh and
# should not be run directly.

. config.h

mkdir -p ${WORK_LW_ORDER_DIR}

GAS_LIST="composite h2o_median o3_median co2_present ch4_present n2o_present o2n2_constant"
# In the microwave we consider limited number of gases
#GAS_LIST="o2n2_constant h2o_median o3_median"

# Loop through the median/present concentrations of each gas and
# reorder
for GAS_SCENARIO in $GAS_LIST
do
    GAS=$(echo ${GAS_SCENARIO} | awk -F_ '{print $1}')
    GAS1=$(echo ${GAS} | awk -F- '{print $1}')

    if [ "$GAS1" = composite ]
    then
	INPUT=${WELL_MIXED_LW_SPECTRA}
    elif [ "$GAS1" = o2n2 ]
    then
	INPUT=${WELL_MIXED_LW_SPECTRA_O2N2}
    else
	INPUT=${MMM_LW_SPECTRA_DIR}/ckdmip_mmm_lw_spectra_${GAS_SCENARIO}.h5
    fi

    for BANDSTRUCT in ${BAND_STRUCTURE}
    do

	OUTPUT=${WORK_LW_ORDER_DIR}/lw_order_${BANDSTRUCT}_${GAS}.h5

	if [ ! -f ${OUTPUT} ]
	then
	    ${BANNER} Reordering ${GAS}, band structure ${BANDSTRUCT}
	    if [ "$BANDSTRUCT" = narrow ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    "wavenumber1=$WN1_LW_NARROW" "wavenumber2=$WN2_LW_NARROW"
	    elif [ "$BANDSTRUCT" = wide ]
	    then
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    "wavenumber1=$WN1_LW_WIDE" "wavenumber2=$WN2_LW_WIDE"
	    elif [ "$BANDSTRUCT" = fsck ]
	    then
		# Assuming FSCK
		${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT
	    else
		if [ "$WN1_LW_CUSTOM" ]
		then
		    ${REORDER_SPECTRUM} iprofile=0 input=$INPUT output=$OUTPUT \
		    ${OPTIONS} \
		    "wavenumber1=$WN1_LW_CUSTOM" "wavenumber2=$WN2_LW_CUSTOM"
		else
		    ${BANNER_ERROR} "Band structure \"$BANDSTRUCT\" not understood"
		    exit 1
		fi
	    fi
	else
	    ${BANNER_SKIP} Skipping reordering of ${GAS} as file present: ${OUTPUT}
	fi

    done
done
