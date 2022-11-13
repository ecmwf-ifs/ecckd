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
# This script finds longwave g points. It is called from do_all_lw.sh
# and should not be run directly. It requires variables TOLERANCE to
# contain a list of heating-rate tolerances (in K d-1) for each g
# point to be used, APPLICATION to be ONE of "climate", "global-nwp"
# or "limited-area-nwp", and BAND_STRUCTURE to contain a list from
# "fsck", "wide" and "narrow". These variables may be provided either
# as environment variables, or should be set in a script provided as
# the first argument to this script that is sourced at the start.

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

# Optional additional arguments
#EXTRA_ARGS="averaging_method=logarithmic"

# Create output directory, if needed
mkdir -p ${WORK_LW_GPOINTS_DIR}

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do

# Create configuration file
if [ "$APP" = climate0 ]
then
    cat > config_find_g_points_lw_climate0.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for climate applications. Absorption spectra for each
# gas are reordered using present-day concentrations, setting the
# concentrations of all other gases also to the minimum over the range
# to be simulated. All well-mixed gases are combined into a single
# hybrid gas.

append_path "${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}"
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.01
flux_weight 0.0
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases composite h2o o3

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_lw_spectra_h2o_median.h5
  reordering_input lw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_lw_spectra_composite_minimum.h5
            ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end h2o

\begin o3
  input ckdmip_mmm_lw_spectra_o3_median.h5
  reordering_input lw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_lw_spectra_composite_minimum.h5
            ckdmip_mmm_lw_spectra_h2o_minimum.h5"
\end o3

\begin composite
  input ckdmip_mmm_lw_spectra_composite_present.h5
  reordering_input lw_order_${BANDSTRUCT}_composite.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
            ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end composite

EOF
elif [ "$APP" = climate2 ]
then
    cat > config_find_g_points_lw_climate2.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for climate applications. Absorption spectra for each
# gas are reordered using present-day concentrations, setting the
# concentrations of all other gases also to the minimum over the range
# to be simulated. 

append_path "${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}"
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.01 
flux_weight 0.0
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases ch4 n2o co2 h2o o3

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_lw_spectra_h2o_median.h5
  reordering_input lw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_ch4_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5"
  background_conc -1 180e-6 350e-9 190e-9
\end h2o

\begin o3
  input ckdmip_mmm_lw_spectra_o3_median.h5
  reordering_input lw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_ch4_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5"
  background_conc -1 180e-6 350e-9 190e-9
\end o3

\begin co2
  input ckdmip_mmm_lw_spectra_co2_present.h5
  reordering_input lw_order_${BANDSTRUCT}_co2.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_ch4_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5"
  background_conc -1 -1 350e-9 190e-9
\end co2

\begin ch4
  input ckdmip_mmm_lw_spectra_ch4_present.h5
  reordering_input lw_order_${BANDSTRUCT}_ch4.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5"
  background_conc -1 -1 180e-6 190e-9
\end ch4

\begin n2o
  input ckdmip_mmm_lw_spectra_n2o_present.h5
  reordering_input lw_order_${BANDSTRUCT}_n2o.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_ch4_present.h5"
  background_conc -1 -1 180e-6 350e-9
\end n2o

EOF
elif [ "$APP" = climate ]
then
    cat > config_find_g_points_lw_climate.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for climate applications. Absorption spectra for each
# gas are reordered using present-day concentrations, setting the
# concentrations of all other gases also to the minimum over the range
# to be simulated. 

append_path "${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}"
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.01 
flux_weight 0.0
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases h2o o3 ch4 n2o co2

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_lw_spectra_h2o_median.h5
  reordering_input lw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_lw_spectra_composite_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end h2o

\begin o3
  input ckdmip_mmm_lw_spectra_o3_median.h5
  reordering_input lw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_lw_spectra_composite_minimum.h5
ckdmip_mmm_lw_spectra_h2o_minimum.h5"
\end o3

\begin co2
  input ckdmip_mmm_lw_spectra_co2_present.h5
  reordering_input lw_order_${BANDSTRUCT}_co2.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_ch4_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5
ckdmip_mmm_lw_spectra_o2_constant.h5
ckdmip_mmm_lw_spectra_n2_constant.h5"
  background_conc -1 -1 350e-9 190e-9 -1 -1
\end co2

\begin ch4
  input ckdmip_mmm_lw_spectra_ch4_present.h5
  reordering_input lw_order_${BANDSTRUCT}_ch4.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5
ckdmip_mmm_lw_spectra_o2_constant.h5
ckdmip_mmm_lw_spectra_n2_constant.h5"
  background_conc -1 -1 180e-6 190e-9 -1 -1
\end ch4

\begin n2o
  input ckdmip_mmm_lw_spectra_n2o_present.h5
  reordering_input lw_order_${BANDSTRUCT}_n2o.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_ch4_present.h5
ckdmip_mmm_lw_spectra_o2_constant.h5
ckdmip_mmm_lw_spectra_n2_constant.h5"
  background_conc -1 -1 180e-6 350e-9 -1 -1
\end n2o

\begin o2n2
  input ckdmip_mmm_lw_spectra_o2n2_constant.h5
  reordering_input lw_order_${BANDSTRUCT}_o2n2.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5
ckdmip_mmm_lw_spectra_co2_present.h5
ckdmip_mmm_lw_spectra_ch4_present.h5
ckdmip_mmm_lw_spectra_n2o_present.h5"
  background_conc -1 -1 180e-6 350e-9 190e-9
\end o2n2

EOF
elif [ "$APP" = nwp ]
then
    cat > config_find_g_points_lw_nwp.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for NWP around the year 2020. Absorption spectra for
# each gas are reordered using present-day concentrations, setting the
# concentrations of all well-mixed gases also to present-day levels. All
# well-mixed gases ar combined into a single hybrid gas.

append_path "${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}"
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.015
flux_weight 0.0
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases composite h2o o3

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_lw_spectra_h2o_median.h5
  reordering_input lw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_lw_spectra_composite_present.h5
            ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end h2o

\begin o3
  input ckdmip_mmm_lw_spectra_o3_median.h5
  reordering_input lw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_lw_spectra_composite_present.h5
            ckdmip_mmm_lw_spectra_h2o_minimum.h5"
\end o3

\begin composite
  input ckdmip_mmm_lw_spectra_composite_present.h5
  reordering_input lw_order_${BANDSTRUCT}_composite.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
            ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end composite

EOF
elif [ "$APP" = nwp-microwave ]
then
    cat > config_find_g_points_lw_nwp-microwave.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for NWP around the year 2020. Absorption spectra for
# each gas are reordered using present-day concentrations, setting the
# concentrations of all well-mixed gases also to present-day levels. All
# well-mixed gases ar combined into a single hybrid gas.

append_path "${MMM_LW_SPECTRA_DIR}:${WORK_LW_SPECTRA_DIR}:${WORK_LW_ORDER_DIR}"
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.015
flux_weight 0.0
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases o2n2 h2o o3

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_lw_spectra_h2o_median.h5
  reordering_input lw_order_${BANDSTRUCT}_h2o.h5
  background_input "ckdmip_mmm_lw_spectra_o2n2_constant.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end h2o

\begin o2n2
  input "ckdmip_mmm_lw_spectra_o2_constant.h5
ckdmip_mmm_lw_spectra_n2_constant.h5"
  reordering_input lw_order_${BANDSTRUCT}_o2n2.h5
  background_input "ckdmip_mmm_lw_spectra_h2o_minimum.h5
ckdmip_mmm_lw_spectra_o3_minimum.h5"
\end o2n2

\begin o3
  input ckdmip_mmm_lw_spectra_o3_median.h5
  reordering_input lw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_lw_spectra_o2n2_constant.h5
ckdmip_mmm_lw_spectra_h2o_minimum.h5"
\end o3

EOF
else
    ${BANNER_ERROR} 'APP "'$APP'" not understood'
    exit 1
fi

    for TOL in $TOLERANCE
    do

	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}

	${BANNER} Finding g-points: $MODEL_CODE

	# When the number of g points is larger than around 30, the
	# FSCK band structure requires at least 3 methane-specific g
	# points or the methane forcing is not very accurate
	CH4_MIN_G_POINTS=
	if [ "$BANDSTRUCT" = fsck -a $(echo "$TOL < 0.018" | bc -l) = 1 ]
	then
	    CH4_MIN_G_POINTS="ch4.min_g_points=3"
	fi

	# Split the water vapour base g point if number of g points
	# larger than around 22
	H2O_BASE_SPLIT=
	if [ "$BANDSTRUCT" = fsck -a $(echo "$TOL < 0.035" | bc -l) = 1 ]
	then
	    H2O_BASE_SPLIT="h2o.base_split=2"
	fi

	${FIND_G_POINTS} \
	    heating_rate_tolerance=${TOL} \
	    output=${WORK_LW_GPOINTS_DIR}/${ECCKD_PREFIX}_lw_gpoints_${MODEL_CODE}.h5 \
	    $EXTRA_ARGS $H2O_BASE_SPLIT $CH4_MIN_G_POINTS config_find_g_points_lw_${APP}.cfg \
	    |& tee ${WORK_LW_GPOINTS_DIR}/${ECCKD_PREFIX}_lw_gpoints_${MODEL_CODE}.log
	test "${PIPESTATUS[0]}" -eq 0
    done
done
