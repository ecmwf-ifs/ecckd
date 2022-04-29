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
# This script optimizes CKD shortwave look-up tables. The input
# requirements are the same as find_g_points_sw.sh. It is called from
# do_all_sw.sh and should not be run directly.

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

if [ -z "$1" ]
then
    OPTIMIZE_MODE_LIST=$APP
else
    OPTIMIZE_MODE_LIST=$@
fi

for OPTIMIZE_MODE in $OPTIMIZE_MODE_LIST
do

# ECCKD 0.6
COMMON_OPTIONS="prior_error=8.0 broadband_weight=0.5 flux_weight=0.1 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=1500"

# ECCKD 0.7
# Too much broadband weight in GB, because the error is dominated by the first band
#COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.3 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000 spectral_boundary_weight=0.1"
COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.3 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000 spectral_boundary_weight=0.0"
# Tried with "double" and "rgb"
COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.3 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000"

COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.4 flux_profile_weight=0.1 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000"
# Tried max_no_rayleigh_wavenumber=15000 but TOA fluxes were improved
# at the expense of surface fluxes, and it was a disaster for CH4 and
# N2O

# Testing values of spectral_boundary_weight of 0, 0.05 and 0.1
# clearly show that 0 is best for broadband fluxes because the
# permitted trade of errors between g points. The nonzero values lead
# to slightly better fluxes per g point, but for cloud radiative
# effect there is nothing to distinguish the three values.

# Test for ecCKD-1.2/fine - slight improvement with the first but the
# second mucks up CH4/N2O forcing - need to have adaptive errors
# according to the range of absorptions being covered by a single g
# point:
#COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.4 flux_weight=0.5 flux_profile_weight=0.2 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000 spectral_boundary_weight=0.05"
#COMMON_OPTIONS="prior_error=2.0 broadband_weight=0.2 flux_weight=0.5 flux_profile_weight=0.2 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 max_iterations=2000 spectral_boundary_weight=0.05"

case "$OPTIMIZE_MODE" in

    nwp)
	# NWP optimization is a single step using present-day gas
	# concentrations
	TRAINING=ckdmip_evaluation1_sw_fluxes_present.h5
	GASLIST="h2o o3 composite"
	INDIR=${WORK_SW_RAW_CKD_DIR}
	OUTDIR=${WORK_SW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=ckd-definition
        SPECIFIC_OPTIONS="convergence_criterion=0.01"
	;;

    climate-base)
	# In the first "base" pass, the properties of composite, H2O,
	# O3 and CO2 are optimized. Then the minor
	# gases in the following passes.

	# First pass
	TRAINING="ckdmip_evaluation1_sw_fluxes_co2-180.h5  ckdmip_evaluation1_sw_fluxes_co2-280.h5
                  ckdmip_evaluation1_sw_fluxes_present.h5  ckdmip_evaluation1_sw_fluxes_co2-560.h5
                  ckdmip_evaluation1_sw_fluxes_co2-1120.h5 ckdmip_evaluation1_sw_fluxes_co2-2240.h5"
	GASLIST="o2n2 h2o o3 co2"
	INDIR=${WORK_SW_RAW_CKD_DIR}
	OUTDIR=${WORK_SW_RAW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=raw2-ckd-definition
        SPECIFIC_OPTIONS="convergence_criterion=0.01 spectral_boundary_weight=0.0"
	;;

    relative-base)
	# "Relative" treats CH4 and N2O using the relative-linear
	# concentration dependence, using present-day concentrations
	# as reference.  The composite gas in this approach contains
	# constant present-day concentrations of O2, N2, CH4 and N2O.

	# In the first "base" pass, the properties of composite, H2O,
	# O3 and CO2 are optimized, where the "rel" reference datasets
	# have CH4 and N2O set to constant profiles. Then the minor
	# gases in the following passes.

	# First pass
	TRAINING="ckdmip_evaluation1_sw_fluxes_rel-180.h5  ckdmip_evaluation1_sw_fluxes_rel-280.h5
                  ckdmip_evaluation1_sw_fluxes_rel-415.h5  ckdmip_evaluation1_sw_fluxes_rel-560.h5
                  ckdmip_evaluation1_sw_fluxes_rel-1120.h5 ckdmip_evaluation1_sw_fluxes_rel-2240.h5"
	# Optionally add one present-day-like scenario from
	# Evaluation-2 to increase the scope of the training data
	if [ "$TRAINING_BOTH" = yes ]
	then
	    TRAINING="$TRAINING ckdmip_evaluation2_sw_fluxes_rel-415.h5"
	fi
	GASLIST="composite h2o o3 co2"
	INDIR=${WORK_SW_RAW_CKD_DIR}
	OUTDIR=${WORK_SW_RAW_CKD_DIR}
#	INCODE=raw-ckd-definition
	INCODE=scaled-ckd-definition
	OUTCODE=raw2-ckd-definition
# New 21 May
	#OPTIONS="prior_error=8.0 broadband_weight=0.5 flux_weight=0.075 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.01"
        #SPECIFIC_OPTIONS="convergence_criterion=0.01 spectral_boundary_weight=0.05"
        SPECIFIC_OPTIONS="convergence_criterion=0.01 spectral_boundary_weight=0.0"
	#"rayleigh_prior_error=1.0"
	;;

    relative-ch4)
	# Second pass
	TRAINING="ckdmip_evaluation1_sw_fluxes_present.h5  ckdmip_evaluation1_sw_fluxes_ch4-350.h5
                  ckdmip_evaluation1_sw_fluxes_ch4-700.h5  ckdmip_evaluation1_sw_fluxes_ch4-1200.h5
                  ckdmip_evaluation1_sw_fluxes_ch4-2600.h5 ckdmip_evaluation1_sw_fluxes_ch4-3500.h5"
	GASLIST="ch4"
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_rel-415.h5"
	INDIR=${WORK_SW_RAW_CKD_DIR}
	OUTDIR=${WORK_SW_RAW_CKD_DIR}
	INCODE=raw2-ckd-definition
	OUTCODE=raw3-ckd-definition
	#OPTIONS="prior_error=8.0 broadband_weight=0.2 flux_weight=0.02 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.0005"
	#OPTIONS="prior_error=8.0 broadband_weight=0.2 flux_weight=0.03 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.0005"
	#OPTIONS="prior_error=8.0 broadband_weight=0.5 flux_weight=0.03 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.0005"
        SPECIFIC_OPTIONS="convergence_criterion=0.0005"
	;;

    relative-n2o)
	# Third pass
        TRAINING="ckdmip_evaluation1_sw_fluxes_present.h5
                  ckdmip_evaluation1_sw_fluxes_n2o-190.h5    ckdmip_evaluation1_sw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_sw_fluxes_n2o-405.h5    ckdmip_evaluation1_sw_fluxes_n2o-540.h5"
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_rel-415.h5"
	GASLIST="n2o"
	INDIR=${WORK_SW_RAW_CKD_DIR}
	OUTDIR=${WORK_SW_CKD_DIR}
	INCODE=raw3-ckd-definition
	OUTCODE=ckd-definition
	#OPTIONS="prior_error=8.0 broadband_weight=0.2 flux_weight=0.02 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.0005"
	#OPTIONS="prior_error=8.0 broadband_weight=0.2 flux_weight=0.03 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.0005"
	#OPTIONS="prior_error=8.0 broadband_weight=0.5 flux_weight=0.03 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 convergence_criterion=0.0005"
        SPECIFIC_OPTIONS="convergence_criterion=0.0005"
	;;

    relative-minor)
	# Second pass doing both CH4 and N2O
	TRAINING="ckdmip_evaluation1_sw_fluxes_present.h5  ckdmip_evaluation1_sw_fluxes_ch4-350.h5
                  ckdmip_evaluation1_sw_fluxes_ch4-700.h5  ckdmip_evaluation1_sw_fluxes_ch4-1200.h5
                  ckdmip_evaluation1_sw_fluxes_ch4-2600.h5 ckdmip_evaluation1_sw_fluxes_ch4-3500.h5
                  ckdmip_evaluation1_sw_fluxes_n2o-190.h5  ckdmip_evaluation1_sw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_sw_fluxes_n2o-405.h5  ckdmip_evaluation1_sw_fluxes_n2o-540.h5"
	GASLIST="ch4 n2o"
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_sw_fluxes_rel-415.h5"
	INDIR=${WORK_SW_RAW_CKD_DIR}
	OUTDIR=${WORK_SW_CKD_DIR}
	INCODE=raw2-ckd-definition
	OUTCODE=ckd-definition
        SPECIFIC_OPTIONS="convergence_criterion=0.0005"
	;;

    *)
	echo 'Optimize mode "'$OPTIMIZE_MODE'" not understood'
	exit 1
esac


mkdir -p ${OUTDIR}

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do
    if [ "$BANDSTRUCT" = wide ]
    then
	# Map from narrow to wide bands
	BANDMAPPING="band_mapping=0 0 0 1 1 1 1 2 2 3 3 4 4"
    elif [ "$BANDSTRUCT" = fsck ]
    then
	BANDMAPPING="band_mapping=0 0 0 0 0 0 0 0 0 0 0 0 0"
    elif [ "$BANDSTRUCT" = double ]
    then
	#BANDMAPPING="band_mapping=0 0 0 0 0 0 0 0 1 1 1 1 1"
	BANDMAPPING="band_mapping=0 0 0 0 0 0 0 0 0 1 1 1 1"
    elif [ "$BANDSTRUCT" = rgb ]
    then
	# Assume we are training from the sw_fluxes-rgb LBL files
	BANDMAPPING="band_mapping=0 0 0 0 1 2 3 4 4"
	# Modify training files and directory
	TRAINING_SW_FLUXES_DIR=$(echo $TRAINING_SW_FLUXES_DIR | sed 's|sw_fluxes$|sw_fluxes-rgb|g')
	TRAINING=$(echo $TRAINING | sed 's/sw_fluxes_/sw_fluxes-rgb_/g')
	EXTRA_ARGS=$(echo $EXTRA_ARGS | sed 's/sw_fluxes_/sw_fluxes-rgb_/g')
    elif [ "$BANDSTRUCT" = gb ]
    then
	# Assume we are training from the sw_fluxes-rgb LBL files
	BANDMAPPING="band_mapping=0 0 0 1 1 2 3 4 4"
	# Modify training files and directory
	TRAINING_SW_FLUXES_DIR=$(echo $TRAINING_SW_FLUXES_DIR | sed 's|sw_fluxes$|sw_fluxes-rgb|g')
	TRAINING=$(echo $TRAINING | sed 's/sw_fluxes_/sw_fluxes-rgb_/g')
	EXTRA_ARGS=$(echo $EXTRA_ARGS | sed 's/sw_fluxes_/sw_fluxes-rgb_/g')
    elif [ "$BANDSTRUCT" = fine ]
    then
	# Assume we are training from the sw_fluxes-fine LBL files
	BANDMAPPING="band_mapping=0 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24"
	# Modify training files and directory
	TRAINING_SW_FLUXES_DIR=$(echo $TRAINING_SW_FLUXES_DIR | sed 's|sw_fluxes$|sw_fluxes-fine|g')
	TRAINING=$(echo $TRAINING | sed 's/sw_fluxes_/sw_fluxes-fine_/g')
	EXTRA_ARGS=$(echo $EXTRA_ARGS | sed 's/sw_fluxes_/sw_fluxes-fine_/g')
    elif [ "$BANDSTRUCT" = window ]
    then
	# Assume we are training from the sw_fluxes-fine LBL files
	BANDMAPPING="band_mapping=0 0 1 2 3 4 5 5 5 6 6 7 7 8 8 9 10 11 12 13 14 15 16 17 18 18"
	# Modify training files and directory
	TRAINING_SW_FLUXES_DIR=$(echo $TRAINING_SW_FLUXES_DIR | sed 's|sw_fluxes$|sw_fluxes-fine|g')
	TRAINING=$(echo $TRAINING | sed 's/sw_fluxes_/sw_fluxes-fine_/g')
	EXTRA_ARGS=$(echo $EXTRA_ARGS | sed 's/sw_fluxes_/sw_fluxes-fine_/g')
    else
	# narrow
	BANDMAPPING="band_mapping=0 1 2 3 4 5 6 7 8 9 10 11 12"
    fi

    cat > config_optimize_lut_${OPTIMIZE_MODE}.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${TRAINING_SW_FLUXES_DIR}:${WORK_SW_LBL_FLUX_DIR}
gases $GASLIST
training_input "$TRAINING"
$(echo $BANDMAPPING | tr "=" " ")
$(echo $COMMON_OPTIONS $SPECIFIC_OPTIONS | tr "= " " \n")
EOF

    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	${BANNER} Optimizing CKD model: $MODEL_CODE

	GPOINTFILE=${WORK_SW_GPOINTS_DIR}/${ECCKD_PREFIX}_sw_gpoints_${MODEL_CODE}.h5
	INPUT=${INDIR}/${ECCKD_PREFIX}_sw_${INCODE}_${MODEL_CODE}.nc
	OUTPUT=${OUTDIR}/${ECCKD_PREFIX}_sw_${OUTCODE}_${MODEL_CODE}.nc
	LOG=${OUTDIR}/${ECCKD_PREFIX}_sw_${OUTCODE}_${MODEL_CODE}.log

	$OPTIMIZE_LUT \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    gpointfile=${GPOINTFILE} \
	    model_id=sw_${APPLICATION}_${BANDSTRUCT}-tol${TOL} \
	    $EXTRA_ARGS \
	    config_optimize_lut_${OPTIMIZE_MODE}.cfg \
	    |& tee $LOG
	test "${PIPESTATUS[0]}" -eq 0
    done
done

done
