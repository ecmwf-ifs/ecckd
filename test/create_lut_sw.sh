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
# This script creates raw shortwave CKD look-up tables. It is called
# from do_all_sw.sh and should not be run directly. The input
# requirements are the same as find_g_points_sw.sh

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

# Create configuration file
if [ "$APP" = nwp ]
then
    cat > config_create_lut_nwp.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_SW_GPOINTS_DIR}:${IDEALIZED_SW_SPECTRA_DIR}:${MMM_CONC_DIR}
gases composite h2o o3
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_sw_spectra_h2o_constant-a.h5
ckdmip_idealized_sw_spectra_h2o_constant-b.h5
ckdmip_idealized_sw_spectra_h2o_constant-c.h5
ckdmip_idealized_sw_spectra_h2o_constant-d.h5
ckdmip_idealized_sw_spectra_h2o_constant-e.h5
ckdmip_idealized_sw_spectra_h2o_constant-f.h5
ckdmip_idealized_sw_spectra_h2o_constant-g.h5
ckdmip_idealized_sw_spectra_h2o_constant-h.h5
ckdmip_idealized_sw_spectra_h2o_constant-i.h5
ckdmip_idealized_sw_spectra_h2o_constant-j.h5
ckdmip_idealized_sw_spectra_h2o_constant-k.h5
ckdmip_idealized_sw_spectra_h2o_constant-l.h5"
\end h2o

\begin o3
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_o3_constant.h5
\end o3

\begin composite
  conc_dependence none
  input "ckdmip_idealized_sw_spectra_co2_constant.h5
ckdmip_idealized_sw_spectra_ch4_constant.h5
ckdmip_idealized_sw_spectra_n2o_constant.h5
ckdmip_idealized_sw_spectra_o2_constant.h5
ckdmip_idealized_sw_spectra_n2_constant.h5"
  conc_input ckdmip_mmm_concentrations.nc
  iprofile 0
\end composite

EOF
elif [ "$APP" = climate-linear ]
then
    cat > config_create_lut_climate-linear.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_SW_GPOINTS_DIR}:${IDEALIZED_SW_SPECTRA_DIR}:${MMM_CONC_DIR}
gases o2n2 h2o o3 co2 ch4 n2o
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_sw_spectra_h2o_constant-a.h5
ckdmip_idealized_sw_spectra_h2o_constant-b.h5
ckdmip_idealized_sw_spectra_h2o_constant-c.h5
ckdmip_idealized_sw_spectra_h2o_constant-d.h5
ckdmip_idealized_sw_spectra_h2o_constant-e.h5
ckdmip_idealized_sw_spectra_h2o_constant-f.h5
ckdmip_idealized_sw_spectra_h2o_constant-g.h5
ckdmip_idealized_sw_spectra_h2o_constant-h.h5
ckdmip_idealized_sw_spectra_h2o_constant-i.h5
ckdmip_idealized_sw_spectra_h2o_constant-j.h5
ckdmip_idealized_sw_spectra_h2o_constant-k.h5
ckdmip_idealized_sw_spectra_h2o_constant-l.h5"
\end h2o

\begin o3
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_o3_constant.h5
\end o3

\begin co2
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_co2_constant.h5
\end co2

\begin ch4
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_ch4_constant.h5
\end ch4

\begin n2o
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_n2o_constant.h5
\end n2o

\begin o2n2
  conc_dependence none
  input "ckdmip_idealized_sw_spectra_o2_constant.h5
ckdmip_idealized_sw_spectra_n2_constant.h5"
  conc_input ckdmip_mmm_concentrations.nc
  iprofile 0
\end o2n2
EOF

elif [ "$APP" = climate ]
then
    cat > config_create_lut_climate.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${WORK_SW_GPOINTS_DIR}:${IDEALIZED_SW_SPECTRA_DIR}:${MMM_CONC_DIR}
gases composite h2o o3 co2 ch4 n2o
#base_wavenumber_boundary 7150 10500
#base_wavenumber_boundary 4000 6000 8000 10000 12000
\begin h2o
  conc_dependence lut
  input "ckdmip_idealized_sw_spectra_h2o_constant-a.h5
ckdmip_idealized_sw_spectra_h2o_constant-b.h5
ckdmip_idealized_sw_spectra_h2o_constant-c.h5
ckdmip_idealized_sw_spectra_h2o_constant-d.h5
ckdmip_idealized_sw_spectra_h2o_constant-e.h5
ckdmip_idealized_sw_spectra_h2o_constant-f.h5
ckdmip_idealized_sw_spectra_h2o_constant-g.h5
ckdmip_idealized_sw_spectra_h2o_constant-h.h5
ckdmip_idealized_sw_spectra_h2o_constant-i.h5
ckdmip_idealized_sw_spectra_h2o_constant-j.h5
ckdmip_idealized_sw_spectra_h2o_constant-k.h5
ckdmip_idealized_sw_spectra_h2o_constant-l.h5"
\end h2o

\begin o3
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_o3_constant.h5
\end o3

\begin co2
  conc_dependence linear
  input ckdmip_idealized_sw_spectra_co2_constant.h5
\end co2

\begin ch4
  conc_dependence relative-linear
  input ckdmip_idealized_sw_spectra_ch4_constant.h5
  reference_conc 1921e-9
\end ch4

\begin n2o
  conc_dependence relative-linear
  input ckdmip_idealized_sw_spectra_n2o_constant.h5
  reference_conc 332e-9
\end n2o

\begin composite
  conc_dependence none
  input "ckdmip_idealized_sw_spectra_o2_constant.h5
ckdmip_idealized_sw_spectra_n2_constant.h5
ckdmip_idealized_sw_spectra_n2o_constant.h5
ckdmip_idealized_sw_spectra_ch4_constant.h5"
  conc_input ckdmip_mmm-const_concentrations.nc
  iprofile 0
\end composite
EOF

else
    ${BANNER_ERROR} 'APP "'$APP'" not understood'
    exit 1
fi

mkdir -p ${WORK_SW_RAW_CKD_DIR}

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do

    if [ "$BANDSTRUCT" = gb ]
    then
	# Split the base spectral interval for the 8000-16650 band at
	# the red/infrared boundary of 700 nm
	EXTRA_ARGS="base_wavenumber_boundary=14300"
    else
	EXTRA_ARGS=
    fi

    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	${BANNER} Creating raw CKD model: $MODEL_CODE

	INPUT=${ECCKD_PREFIX}_sw_gpoints_${MODEL_CODE}.h5
	OUTPUT=${WORK_SW_RAW_CKD_DIR}/${ECCKD_PREFIX}_sw_raw-ckd-definition_${MODEL_CODE}.nc
	$CREATE_LOOK_UP_TABLE \
	    ssi=$MMM_SW_SSI \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    temperature_stride=1 $EXTRA_ARGS \
	    config_create_lut_${APP}.cfg \
	    |& tee ${WORK_SW_RAW_CKD_DIR}/${ECCKD_PREFIX}_sw_raw-ckd-definition_${MODEL_CODE}.log
	test "${PIPESTATUS[0]}" -eq 0
    done
done
