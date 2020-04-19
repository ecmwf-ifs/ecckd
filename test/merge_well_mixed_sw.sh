#!/bin/bash
# Create composite gas containing all gases except H2O and O3.  Note
# that the "--merge-only" option of ckdmip_sw doesn't seem to work
# properly so we use ckdmip_lw.

. config.h

mkdir -p ${WORK_SW_SPECTRA_DIR}

# Create merge for NWP application
if [ ! -f ${WELL_MIXED_SW_SPECTRA} ]
then
    ${BANNER} Merging spectra of well-mixed gases for present-day conditions

    GAS_FILES=
    for GAS_SCENARIO in o2_constant n2_constant co2_present ch4_present n2o_present cfc11_present-equivalent cfc12_present
    do
	GAS_FILES="${GAS_FILES} ${MMM_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_${GAS_SCENARIO}.h5"
    done
    ${CKDMIP_LW} \
	--merge-only ${GAS_FILES} \
	--output ${WELL_MIXED_SW_SPECTRA}

else
    ${BANNER_SKIP} Skipping merge of well-mixed gases for present-day conditions as file present: ${WELL_MIXED_SW_SPECTRA}
fi

# Create merge for climate applications
if [ ! -f ${WELL_MIXED_SW_SPECTRA_MINIMUM} ]
then
    ${BANNER} Merging spectra of well-mixed gases for climate minimum conditions

    GAS_FILES="${MMM_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_o2_constant.h5
${MMM_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_n2_constant.h5
--conc 180e-6 ${MMM_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_co2_present.h5
--conc 350e-9 ${MMM_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_ch4_present.h5
--conc 190e-9 ${MMM_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_n2o_present.h5"
    ${CKDMIP_LW} \
	--merge-only ${GAS_FILES} \
	--output ${WELL_MIXED_SW_SPECTRA_MINIMUM}
else
    ${BANNER_SKIP} Skipping merge of well-mixed gases for climate minimum as file present: ${WELL_MIXED_SW_SPECTRA_MINIMUM}
fi
