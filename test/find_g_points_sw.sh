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
# This script finds shortwave g points. It is called from do_all_sw.sh
# and should not be run directly.  It requires variables TOLERANCE to
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
EXTRA_ARGS="averaging_method=total-transmission max_no_rayleigh_wavenumber=10000"

# Create output directory, if needed
mkdir -p ${WORK_SW_GPOINTS_DIR}

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do

O3_MIN_G_POINTS=""
#CO2_MIN_G_POINTS=""
H2O_SPLIT=""
CH4_MIN_G_POINTS=""
N2O_MIN_G_POINTS=""

# Note that min_g_points lists the minimum number of g points per
# band, where any further bands are assumed to have a minimum of 1
if [ "$BANDSTRUCT" = "rgb" -o "$BANDSTRUCT" = "gb" ]
then
    # Need at least 3 g-points for ozone in the UV band
    O3_MIN_G_POINTS="min_g_points 1 1 1 1 3"
    # CO2 minimum g points is better specified as a function of TOL later
    #CO2_MIN_G_POINTS="min_g_points 4 1 1 1 1"
elif [ "$BANDSTRUCT" = "fine" ]
then
    CH4_MIN_G_POINTS="min_g_points 2"
    N2O_MIN_G_POINTS="min_g_points 3"
    O3_MIN_G_POINTS="min_g_points 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 4"
elif [ "$BANDSTRUCT" = "vfine" ]
then
    CH4_MIN_G_POINTS="min_g_points 2"
    N2O_MIN_G_POINTS="min_g_points 3"
    O3_MIN_G_POINTS="min_g_points 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5"
elif [ "$BANDSTRUCT" = "window" ]
then
    CH4_MIN_G_POINTS="min_g_points 2"
    N2O_MIN_G_POINTS="min_g_points 2"
    #O3_MIN_G_POINTS=""
    O3_MIN_G_POINTS="min_g_points 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 4"
fi


# Create configuration file
if [ "$APP" = climate0 ]
then
    cat > config_find_g_points_sw_climate.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for climate applications. Absorption spectra for each
# gas are reordered using present-day concentrations, setting the
# concentrations of all other gases also to the minimum over the range
# to be simulated. All well-mixed gases are combined into a single
# hybrid gas.

append_path "${MMM_SW_SPECTRA_DIR}:${WORK_SW_SPECTRA_DIR}:${WORK_SW_ORDER_DIR}"
ssi $MMM_SW_SSI
iprofile 0
averaging_method "total-transmission"
tolerance_tolerance 0.01 
#flux_weight 0.02
flux_weight 0.0002
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases composite h2o o3

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_median.h5
  reordering_input sw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_sw_spectra_composite_minimum.h5
            ckdmip_mmm_sw_spectra_o3_minimum.h5"
  min_scaling 0.1
  max_scaling 10.0
\end h2o

\begin o3
  input ckdmip_mmm_sw_spectra_o3_median.h5
  reordering_input sw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_sw_spectra_composite_minimum.h5
            ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5"
  min_scaling 0.5
  max_scaling 2.0
  $O3_MIN_G_POINTS
\end o3

\begin composite
  input ckdmip_mmm_sw_spectra_composite_present.h5
  reordering_input sw_order_${BANDSTRUCT}_composite.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
            ckdmip_mmm_sw_spectra_o3_minimum.h5"
\end composite

EOF

elif [ "$APP" = climate ]
then
    cat > config_find_g_points_sw_climate.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for climate applications. Absorption spectra for each
# gas are reordered using present-day concentrations, setting the
# concentrations of all other gases also to the minimum over the range
# to be simulated. 

append_path "${MMM_SW_SPECTRA_DIR}:${WORK_SW_SPECTRA_DIR}:${WORK_SW_ORDER_DIR}"
ssi $MMM_SW_SSI
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.01
flux_weight 0.001 # Previous
#flux_weight 0.0002
#flux_weight 0.2
#flux_weight 0.1
min_pressure ${MIN_PRESSURE}
max_iterations 60

gases h2o o3 ch4 n2o co2 o2n2

#cloud liquidcloud
\begin liquidcloud
  reordering_input sw_order_${BANDSTRUCT}_cloud.h5
  max_reflectance_range 0.34
\end liquidcloud

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_median.h5
  reordering_input sw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_sw_spectra_composite_minimum.h5
ckdmip_mmm_sw_spectra_o3_minimum.h5"
  min_scaling 0.1
  max_scaling 10.0
\end h2o

\begin o3
  input ckdmip_mmm_sw_spectra_o3_median.h5
  reordering_input sw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_sw_spectra_composite_minimum.h5
ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5"
  min_scaling 0.5
  max_scaling 2.0
  $O3_MIN_G_POINTS
\end o3

\begin co2
  input ckdmip_mmm_sw_spectra_co2_present.h5
  reordering_input sw_order_${BANDSTRUCT}_co2.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
ckdmip_mmm_sw_spectra_o3_minimum.h5
ckdmip_mmm_sw_spectra_ch4_present.h5
ckdmip_mmm_sw_spectra_n2o_present.h5
ckdmip_mmm_sw_spectra_o2n2_constant.h5"
  background_conc -1 -1 350e-9 190e-9 -1
  min_scaling 0.43
  max_scaling 5.4
  $CO2_MIN_G_POINTS
\end co2

\begin ch4
  input ckdmip_mmm_sw_spectra_ch4_present.h5
  reordering_input sw_order_${BANDSTRUCT}_ch4.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
ckdmip_mmm_sw_spectra_o3_minimum.h5
ckdmip_mmm_sw_spectra_co2_present.h5
ckdmip_mmm_sw_spectra_n2o_present.h5
ckdmip_mmm_sw_spectra_o2n2_constant.h5"
  background_conc -1 -1 180e-6 190e-9 -1
  min_scaling 1.8
  max_scaling 0.18
  $CH4_MIN_G_POINTS
\end ch4

\begin n2o
  input ckdmip_mmm_sw_spectra_n2o_present.h5
  reordering_input sw_order_${BANDSTRUCT}_n2o.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
ckdmip_mmm_sw_spectra_o3_minimum.h5
ckdmip_mmm_sw_spectra_co2_present.h5
ckdmip_mmm_sw_spectra_ch4_present.h5
ckdmip_mmm_sw_spectra_o2n2_constant.h5"
  background_conc -1 -1 180e-6 350e-9 -1
  min_scaling 1.6
  max_scaling 0.57
  $N2O_MIN_G_POINTS
\end n2o

\begin o2n2
  input ckdmip_mmm_sw_spectra_o2n2_constant.h5
  reordering_input sw_order_${BANDSTRUCT}_o2n2.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
ckdmip_mmm_sw_spectra_o3_minimum.h5
ckdmip_mmm_sw_spectra_co2_present.h5
ckdmip_mmm_sw_spectra_ch4_present.h5
ckdmip_mmm_sw_spectra_n2o_present.h5"
  background_conc -1 -1 180e-6 350e-9 190e-9
\end o2n2

EOF
elif [ "$APP" = nwp ]
then
    cat > config_find_g_points_sw_nwp.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
# Configuration file for generating g points for a radiation scheme to
# be used only for NWP around the year 2020. Absorption spectra for
# each gas are reordered using present-day concentrations, setting the
# concentrations of all well-mixed gases also to present-day levels. All
# well-mixed gases ar combined into a single hybrid gas.

append_path "${MMM_SW_SPECTRA_DIR}:${WORK_SW_SPECTRA_DIR}:${WORK_SW_ORDER_DIR}"
ssi $MMM_SW_SSI
iprofile 0
averaging_method "transmission"
tolerance_tolerance 0.01
#flux_weight 0.005
#flux_weight 0.001 ### Previous
flux_weight 0.0002
min_pressure ${MIN_PRESSURE}
max_iterations 60

#gases rayleigh composite h2o o3
gases composite h2o o3

\begin h2o
  # Water vapour in median present-day concentrations
  input ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_median.h5
  reordering_input sw_order_${BANDSTRUCT}_h2o.h5
  # Other gases in present-day concentrations, except ozone which uses
  # the minimum concentration
  background_input "ckdmip_mmm_sw_spectra_composite_present.h5
            ckdmip_mmm_sw_spectra_o3_minimum.h5"
  min_scaling 0.1
  max_scaling 10.0
\end h2o

\begin o3
  input ckdmip_mmm_sw_spectra_o3_median.h5

  reordering_input sw_order_${BANDSTRUCT}_o3.h5
  background_input "ckdmip_mmm_sw_spectra_composite_present.h5
            ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5"
  min_scaling 0.5
  max_scaling 2.0
  $O3_MIN_G_POINTS
\end o3

\begin composite
  input ckdmip_mmm_sw_spectra_composite_present.h5
  reordering_input sw_order_${BANDSTRUCT}_composite.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
            ckdmip_mmm_sw_spectra_o3_minimum.h5"
\end composite

\begin rayleigh
  input ckdmip_mmm_sw_spectra_rayleigh_present.h5
  reordering_input sw_order_${BANDSTRUCT}_rayleigh.h5
  background_input "ckdmip_mmm_sw_spectra_h2o${H2OSUFFIX}_minimum.h5
            ckdmip_mmm_sw_spectra_o3_minimum.h5
            ckdmip_mmm_sw_spectra_composite_present.h5"
\end rayleigh

EOF
else
    ${BANNER_ERROR} 'APP "'$APP'" not understood'
    exit 1
fi

    for TOL in $TOLERANCE
    do

	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}

	${BANNER} Finding g-points: $MODEL_CODE

	# Reduce the tolerance for the UV band of the RGB band
	# structure
	TOL_BAND=$TOL
	# if [ "$BANDSTRUCT" = rgb ]
	# then
	#     TOL_BAND="$TOL $TOL $TOL $TOL $(echo $TOL/10.0 | bc -l)"
	# elif [ "$BANDSTRUCT" = gb ]
	# then
	#     TOL_BAND="$TOL $TOL $TOL $TOL $(echo $TOL/15.0 | bc -l)"
	# fi

	# We need at least two N2O-specific g points to get reasonable
	# N2O forcing, but this is not really very important so we
	# only force the extra g point when the total number of g
	# points is larger than around 28
	N2O_MIN_G_POINTS=
	if [ "$BANDSTRUCT" = rgb -a $(echo "$TOL < 0.05" | bc -l) = 1 ]
	then
	    N2O_MIN_G_POINTS="n2o.min_g_points=2"
	fi

	# Optionally split the near-IR water vapour spectrum into 2 or
	# 3 subbands, using bash arrays to allow for correct expansion
	# of h2o.subband_wavenumber_boundary if it contains a space
	unset H2O_SPLIT
	if [ "$BANDSTRUCT" = rgb ]
	then
	    if [ $(echo "$TOL < 0.025" | bc -l) = 1 ]
	    then
		H2O_SPLIT=(h2o.g_split=0.75 "h2o.subband_wavenumber_boundary=3750 5350 7150 8700 10650" co2.min_g_points=7)
	    elif [ $(echo "$TOL < 0.065" | bc -l) = 1 ]
	    then
		H2O_SPLIT=(h2o.g_split=0.7 "h2o.subband_wavenumber_boundary=5350 7150 8700 10650" co2.min_g_points=6 "h2o.max_g_points=256 1")
	    elif [ $(echo "$TOL < 0.15" | bc -l) = 1 ]
	    then
		H2O_SPLIT=(h2o.g_split=0.67 "h2o.subband_wavenumber_boundary=7150 10650" co2.min_g_points=5)
	    elif [ $(echo "$TOL < 0.2" | bc -l) = 1 ]
	    then
		H2O_SPLIT=(h2o.g_split=0.65 h2o.subband_wavenumber_boundary=7150 co2.min_g_points=4)
	    fi
	fi

	# O2N2 composite includes Rayleigh, but is reluctant to split
	# for some reason even when Rayleigh is strong.  This is
	# mainly needed in the visible bands where there is little
	# ozone.
	O2N2_MIN_G_POINTS="o2n2.min_g_points=1"
	if [ "$BANDSTRUCT" = narrow -a $(echo "$TOL < 0.05" | bc -l) = 1 ]
	then
	    O2N2_MIN_G_POINTS="o2n2.min_g_points=1 1 1 1 1 1 1 1 1 2 3 1 1"
	fi

	${FIND_G_POINTS} \
	    "heating_rate_tolerance=${TOL_BAND}" \
	    output=${WORK_SW_GPOINTS_DIR}/${ECCKD_PREFIX}_sw_gpoints_${MODEL_CODE}.h5 \
	    $EXTRA_ARGS $N2O_MIN_G_POINTS "${H2O_SPLIT[@]}" "$O2N2_MIN_G_POINTS" \
	    config_find_g_points_sw_${APP}.cfg \
	    |& tee ${WORK_SW_GPOINTS_DIR}/${ECCKD_PREFIX}_sw_gpoints_${MODEL_CODE}.log
	test "${PIPESTATUS[0]}" -eq 0
    done
done
