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
# This script creates the longwave spectra for composite gases
# containing all gases except H2O and O3.

. config.h

mkdir -p ${WORK_LW_SPECTRA_DIR}

# Create merge for NWP application
if [ ! -f ${WELL_MIXED_LW_SPECTRA} ]
then
    ${BANNER} Merging spectra of well-mixed gases for present-day conditions

    GAS_FILES=
    for GAS_SCENARIO in o2_constant n2_constant co2_present ch4_present n2o_present cfc11_present-equivalent cfc12_present
    do
	GAS_FILES="${GAS_FILES} ${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_${GAS_SCENARIO}.h5"
    done
    ${CKDMIP_LW} \
	--merge-only ${GAS_FILES} \
	--output ${WELL_MIXED_LW_SPECTRA}

else
    ${BANNER_SKIP} Skipping merge of well-mixed gases for present-day conditions as file present: ${WELL_MIXED_LW_SPECTRA}
fi

# Create merge for climate applications
if [ ! -f ${WELL_MIXED_LW_SPECTRA_MINIMUM} ]
then
    ${BANNER} Merging spectra of well-mixed gases for climate minimum conditions

    GAS_FILES="${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_o2_constant.h5
${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_n2_constant.h5
--conc 180e-6 ${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_co2_present.h5
--conc 350e-9 ${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_ch4_present.h5
--conc 190e-9 ${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_n2o_present.h5"
    ${CKDMIP_LW} \
	--merge-only ${GAS_FILES} \
	--output ${WELL_MIXED_LW_SPECTRA_MINIMUM}
else
    ${BANNER_SKIP} Skipping merge of well-mixed gases for climate minimum as file present: ${WELL_MIXED_LW_SPECTRA_MINIMUM}
fi

# Create merge of O2 and N2 for climate applications
if [ ! -f ${WELL_MIXED_LW_SPECTRA_O2N2} ]
then
    ${BANNER} Merging spectra of O2 and N2 for climate minimum conditions

    GAS_FILES="${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_o2_constant.h5
${MMM_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_n2_constant.h5"
    ${CKDMIP_LW} \
	--merge-only ${GAS_FILES} \
	--output ${WELL_MIXED_LW_SPECTRA_O2N2}
else
    ${BANNER_SKIP} Skipping merge of O2 and N2 for climate minimum as file present: ${WELL_MIXED_LW_SPECTRA_O2N2}
fi
