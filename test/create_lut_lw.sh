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
# This script creates raw longwave CKD look-up tables. It is called
# from do_all_lw.sh and should not be run directly. The input
# requirements are the same as find_g_points_lw.sh

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

# Default spectral averaging method matches the transmittance of a
# layer (a tenth of a decade in pressure)
OUTPUT_CODE=_raw

# Alternative is to match the transmittance of three layers (a third
# of a decade in pressure), which seems to be very slightly worse
#EXTRA_ARGS="averaging_method=transmission-3"
# Or a hybrid between logarithmic averaging in the troposphere and
# transmission-3 in the stratosphere, but this appears to totally ruin
# calculations of methane radiative forcing
#EXTRA_ARGS="averaging_method=logarithmic"
#EXTRA_ARGS="averaging_method=hybrid-logarithmic-transmission-3"

# Create configuration file
if [ "$APP" = nwp ]
then
    cat > config_create_lut_nwp.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_LW_GPOINTS_DIR}:${IDEALIZED_LW_SPECTRA_DIR}:${MMM_CONC_DIR}:${WORK_LW_SPECTRA_DIR}
gases composite h2o o3
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-a.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-b.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-c.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-d.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-e.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-f.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-g.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-h.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-i.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-j.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-k.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-l.h5"
\end h2o

\begin o3
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_o3_constant.h5
\end o3

\begin composite
  conc_dependence none
  input "ckdmip_idealized_lw_spectra_co2_constant.h5
ckdmip_idealized_lw_spectra_ch4_constant.h5
ckdmip_idealized_lw_spectra_n2o_constant.h5
ckdmip_idealized_lw_spectra_cfc11_constant-equivalent.h5
ckdmip_idealized_lw_spectra_cfc12_constant.h5
ckdmip_idealized_lw_spectra_o2_constant.h5
ckdmip_idealized_lw_spectra_n2_constant.h5"
  conc_input ckdmip_mmm_concentrations.nc
  iprofile 0
\end composite

EOF
elif [ "$APP" = climate4 ]
then
    cat > config_create_lut_climate.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_LW_GPOINTS_DIR}:${IDEALIZED_LW_SPECTRA_DIR}:${MMM_CONC_DIR}:${WORK_LW_SPECTRA_DIR}
gases o2n2 h2o o3 co2 ch4 n2o cfc11 cfc12
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-a.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-b.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-c.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-d.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-e.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-f.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-g.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-h.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-i.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-j.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-k.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-l.h5"
\end h2o

\begin o3
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_o3_constant.h5
\end o3

\begin co2
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_co2_constant.h5
\end co2

\begin ch4
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_ch4_constant.h5
#  reference_conc 1921e-9
\end ch4

\begin n2o
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_n2o_constant.h5
#  reference_conc 332e-9
\end n2o

\begin cfc11
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_cfc11_constant-equivalent.h5
\end cfc11

\begin cfc12
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_cfc12_constant.h5
\end cfc12

\begin o2n2
  conc_dependence none
  input "ckdmip_idealized_lw_spectra_o2_constant.h5
ckdmip_idealized_lw_spectra_n2_constant.h5"
  conc_input ckdmip_mmm_concentrations.nc
  iprofile 0
\end o2n2
EOF

elif [ "$APP" = climate ]
then
    cat > config_create_lut_climate.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_LW_GPOINTS_DIR}:${IDEALIZED_LW_SPECTRA_DIR}:${MMM_CONC_DIR}:${WORK_LW_SPECTRA_DIR}
gases composite h2o o3 co2 ch4 n2o cfc11 cfc12
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-a.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-b.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-c.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-d.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-e.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-f.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-g.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-h.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-i.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-j.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-k.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-l.h5"
\end h2o

\begin o3
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_o3_constant.h5
\end o3

\begin co2
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_co2_constant.h5
\end co2

\begin ch4
  conc_dependence relative-linear
  input ckdmip_idealized_lw_spectra_ch4_constant.h5
  reference_conc 1921e-9
\end ch4

\begin n2o
  conc_dependence relative-linear
  input ckdmip_idealized_lw_spectra_n2o_constant.h5
  reference_conc 332e-9
\end n2o

\begin cfc11
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_cfc11_constant-equivalent.h5
\end cfc11

\begin cfc12
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_cfc12_constant.h5
\end cfc12

\begin composite
  conc_dependence none
  input "ckdmip_idealized_lw_spectra_o2_constant.h5
ckdmip_idealized_lw_spectra_n2_constant.h5
ckdmip_idealized_lw_spectra_n2o_constant.h5
ckdmip_idealized_lw_spectra_ch4_constant.h5"
  conc_input ckdmip_mmm-const_concentrations.nc
  iprofile 0
\end composite
EOF
elif [ "$APP" = nwp-microwave ]
then
    cat > config_create_lut_nwp-microwave.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_LW_GPOINTS_DIR}:${IDEALIZED_LW_SPECTRA_DIR}:${MMM_CONC_DIR}:${WORK_LW_SPECTRA_DIR}
gases o2n2 h2o o3
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-a.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-b.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-c.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-d.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-e.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-f.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-g.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-h.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-i.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-j.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-k.h5
ckdmip_idealized_lw_spectra_h2o${H2OSUFFIX}_constant-l.h5"
\end h2o

\begin o2n2
  conc_dependence none
  input "ckdmip_idealized_lw_spectra_o2_constant.h5
ckdmip_idealized_lw_spectra_n2_constant.h5"
  conc_input ckdmip_mmm_concentrations.nc
  iprofile 0
\end o2n2

\begin o3
  conc_dependence linear
  input ckdmip_idealized_lw_spectra_o3_constant.h5
\end o3

EOF
else
    ${BANNER_ERROR} 'APP "'$APP'" not understood'
    exit 1
fi

mkdir -p ${WORK_LW_RAW_CKD_DIR}

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	${BANNER} Creating raw CKD model: $MODEL_CODE

	INPUT=${ECCKD_PREFIX}_lw_gpoints_${MODEL_CODE}.h5
	OUTPUT=${WORK_LW_RAW_CKD_DIR}/${ECCKD_PREFIX}_lw${OUTPUT_CODE}-ckd-definition_${MODEL_CODE}.nc
        $CREATE_LOOK_UP_TABLE \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    $EXTRA_ARGS \
	    config_create_lut_${APP}.cfg \
	    |& tee ${WORK_LW_RAW_CKD_DIR}/${ECCKD_PREFIX}_lw${OUTPUT_CODE}-ckd-definition_${MODEL_CODE}.log
	test "${PIPESTATUS[0]}" -eq 0
    done
done
