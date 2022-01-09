#!/bin/bash
# Optimize CKD look-up tables. The input requirements are the same as
# find_g_points_lw.sh

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

#OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_weight=0.05 flux_profile_weight=0.05 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8"

COMMON_OPTIONS="prior_error=8.0 broadband_weight=0.8 flux_profile_weight=0.2 temperature_corr=0.8 pressure_corr=0.8 conc_corr=0.8 spectral_boundary_weight=0.1"

case "$OPTIMIZE_MODE" in

    nwp)
	# NWP optimization is a single step using present-day gas
	# concentrations
	TRAINING=ckdmip_evaluation1_lw_fluxes_present.h5
	GASLIST="h2o o3 composite"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=ckd-definition
	SPECIFIC_OPTIONS="flux_weight=0.2"
	;;

    zero-minor1)
	# Zero-Minor first optimizes O2, N2, H2O, O3 and CO2 assuming
	# the others are zero, then optimizes CH4, N2O and CFC11

	# First pass
	TRAINING="ckdmip_evaluation1_lw_fluxes_5gas-180.h5  ckdmip_evaluation1_lw_fluxes_5gas-280.h5 
                  ckdmip_evaluation1_lw_fluxes_5gas-415.h5  ckdmip_evaluation1_lw_fluxes_5gas-560.h5
                  ckdmip_evaluation1_lw_fluxes_5gas-1120.h5 ckdmip_evaluation1_lw_fluxes_5gas-2240.h5"
	GASLIST="o2n2 h2o o3 co2"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_RAW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=raw2-ckd-definition
	SPECIFIC_OPTIONS="flux_weight=0.2"
	;;

    zero-minor2)
	# Second pass
	TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5  ckdmip_evaluation1_lw_fluxes_ch4-350.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-700.h5  ckdmip_evaluation1_lw_fluxes_ch4-1200.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-2600.h5 ckdmip_evaluation1_lw_fluxes_ch4-3500.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-190.h5  ckdmip_evaluation1_lw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-405.h5  ckdmip_evaluation1_lw_fluxes_n2o-540.h5
                  ckdmip_evaluation1_lw_fluxes_cfc11-0.h5  ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5"
	# Don't bother optimizing CFC12 - more accurate out of the box
	#ckdmip_evaluation1_lw_fluxes_cfc12-0.h5
	#ckdmip_evaluation1_lw_fluxes_cfc12-550.h5
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_5gas-415.h5"
	#EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5"
	GASLIST="ch4 n2o cfc11"
	#    GASLIST="ch4"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_CKD_DIR}
	INCODE=raw2-ckd-definition
	OUTCODE=ckd-definition
	SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.2"
	;;
    
    climate | all-in-one)
	# Optimize all gases together
	GASLIST="composite h2o o3 co2 ch4 n2o cfc11"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=ckd-definition
	TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5  ckdmip_evaluation1_lw_fluxes_co2-180.h5
                  ckdmip_evaluation1_lw_fluxes_co2-280.h5  ckdmip_evaluation1_lw_fluxes_co2-560.h5
                  ckdmip_evaluation1_lw_fluxes_co2-1120.h5 ckdmip_evaluation1_lw_fluxes_co2-2240.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-700.h5  ckdmip_evaluation1_lw_fluxes_ch4-1200.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-2600.h5 ckdmip_evaluation1_lw_fluxes_ch4-3500.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-190.h5  ckdmip_evaluation1_lw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-405.h5  ckdmip_evaluation1_lw_fluxes_n2o-540.h5
                  ckdmip_evaluation1_lw_fluxes_cfc11-0.h5  ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5"
	SPECIFIC_OPTIONS="flux_weight=0.2"
	;;

    all1)
	# Optimize all gases together
	GASLIST="composite h2o o3 co2 ch4 n2o cfc11"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_RAW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=rawall-ckd-definition
	TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5  ckdmip_evaluation1_lw_fluxes_co2-180.h5
                  ckdmip_evaluation1_lw_fluxes_co2-280.h5  ckdmip_evaluation1_lw_fluxes_co2-560.h5
                  ckdmip_evaluation1_lw_fluxes_co2-1120.h5 ckdmip_evaluation1_lw_fluxes_co2-2240.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-700.h5  ckdmip_evaluation1_lw_fluxes_ch4-1200.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-2600.h5 ckdmip_evaluation1_lw_fluxes_ch4-3500.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-190.h5  ckdmip_evaluation1_lw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-405.h5  ckdmip_evaluation1_lw_fluxes_n2o-540.h5
                  ckdmip_evaluation1_lw_fluxes_cfc11-0.h5  ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5"
	SPECIFIC_OPTIONS="convergence_criterion=0.02 flux_weight=0.2"
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
	TRAINING="ckdmip_evaluation1_lw_fluxes_rel-180.h5  ckdmip_evaluation1_lw_fluxes_rel-280.h5
                  ckdmip_evaluation1_lw_fluxes_rel-415.h5  ckdmip_evaluation1_lw_fluxes_rel-560.h5
                  ckdmip_evaluation1_lw_fluxes_rel-1120.h5 ckdmip_evaluation1_lw_fluxes_rel-2240.h5"
	# Optionally add one present-day-like scenario from
	# Evaluation-2 to increase the scope of the training data
	if [ "$TRAINING_BOTH" = yes ]
	then
	    TRAINING="$TRAINING ckdmip_evaluation2_lw_fluxes_rel-415.h5"
	fi

	GASLIST="composite h2o o3 co2"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_RAW_CKD_DIR}
	INCODE=raw-ckd-definition
	OUTCODE=raw2-ckd-definition
	;;

    relative-ch4)
	# Second pass
	TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5  ckdmip_evaluation1_lw_fluxes_ch4-350.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-700.h5  ckdmip_evaluation1_lw_fluxes_ch4-1200.h5
                  ckdmip_evaluation1_lw_fluxes_ch4-2600.h5 ckdmip_evaluation1_lw_fluxes_ch4-3500.h5"
	GASLIST="ch4"
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5"
	INDIR=${WORK_LW_RAW_CKD_DIR}

	OUTDIR=${WORK_LW_RAW_CKD_DIR}
	INCODE=raw2-ckd-definition
	OUTCODE=raw3-ckd-definition
	# Convergence appears to be better in some cases when the
	# penaty for negative optical depths is low (default is 1e4.
	# Note that operationally the optical depth would always be
	# capped to be no lower than 0.
	SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.5 negative_od_penalty=1.0e1"
	#SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.5"
	;;

    relative-minor)
	# Third pass
        TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5  
                  ckdmip_evaluation1_lw_fluxes_n2o-190.h5    ckdmip_evaluation1_lw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-405.h5    ckdmip_evaluation1_lw_fluxes_n2o-540.h5
                  ckdmip_evaluation1_lw_fluxes_cfc11-0.h5    ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5"
	# Don't bother optimizing CFC12 - more accurate out of the box
	#ckdmip_evaluation1_lw_fluxes_cfc12-0.h5
	#ckdmip_evaluation1_lw_fluxes_cfc12-550.h5
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5"
	GASLIST="n2o cfc11"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_CKD_DIR}
	INCODE=raw3-ckd-definition
	OUTCODE=ckd-definition
	SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.5"
	;;

    relative-n2o)
	# Third pass
        TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-190.h5    ckdmip_evaluation1_lw_fluxes_n2o-270.h5
                  ckdmip_evaluation1_lw_fluxes_n2o-405.h5    ckdmip_evaluation1_lw_fluxes_n2o-540.h5"
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5"
	GASLIST="n2o"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_RAW_CKD_DIR}
	INCODE=raw3-ckd-definition
	OUTCODE=raw4-ckd-definition
	SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.5"
	;;

    relative-cfc11)
	# Fourth pass
        TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5
                  ckdmip_evaluation1_lw_fluxes_cfc11-0.h5    ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5"
	# Don't bother optimizing CFC12 - more accurate out of the box
	#ckdmip_evaluation1_lw_fluxes_cfc12-0.h5
	#ckdmip_evaluation1_lw_fluxes_cfc12-550.h5
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5"
	GASLIST="cfc11"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_CKD_DIR}
	INCODE=raw4-ckd-definition
	OUTCODE=ckd-definition
	SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.2"
	;;

    relative-cfc)
	# Fourth pass
        TRAINING="ckdmip_evaluation1_lw_fluxes_present.h5
                  ckdmip_evaluation1_lw_fluxes_cfc11-0.h5    ckdmip_evaluation1_lw_fluxes_cfc11-2000.h5
	          ckdmip_evaluation1_lw_fluxes_cfc12-0.h5    ckdmip_evaluation1_lw_fluxes_cfc12-550.h5"
	EXTRA_ARGS="$EXTRA_ARGS relative_to=ckdmip_evaluation1_lw_fluxes_rel-415.h5"
	GASLIST="cfc11 cfc12"
	INDIR=${WORK_LW_RAW_CKD_DIR}
	OUTDIR=${WORK_LW_CKD_DIR}
	INCODE=raw4-ckd-definition
	OUTCODE=ckd-definition
	SPECIFIC_OPTIONS="convergence_criterion=0.0005 flux_weight=0.2"
	;;

    *)
	echo 'Optimize mode "'$OPTIMIZE_MODE'" not understood'
	exit 1
esac

mkdir -p "${OUTDIR}"

# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do
    if [ "$BANDSTRUCT" = wide ]
    then
	# Map from narrow to wide bands
	BANDMAPPING="band_mapping=0 0 1 1 1 2 2 2 3 3 3 4 4"
    elif [ "$BANDSTRUCT" = fsck ]
    then
	BANDMAPPING="band_mapping=0 0 0 0 0 0 0 0 0 0 0 0 0"
    else
	BANDMAPPING="band_mapping=0 1 2 3 4 5 6 7 8 9 10 11 12"
    fi

    cat > config_optimize_lut_${OPTIMIZE_MODE}.cfg <<EOF
# THIS FILE WAS AUTOMATICALLY GENERATED BY $0
append_path ${TRAINING_LW_FLUXES_DIR}:${WORK_LW_LBL_FLUX_DIR}
gases $GASLIST
training_input "$TRAINING"
$(echo $BANDMAPPING | tr "=" " ")
$(echo $COMMON_OPTIONS $SPECIFIC_OPTIONS | tr "= " " \n")
EOF

    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	${BANNER} Optimizing CKD model: $MODEL_CODE

	GPOINTFILE=${WORK_LW_GPOINTS_DIR}/${ECCKD_PREFIX}_lw_gpoints_${MODEL_CODE}.h5
	INPUT=${INDIR}/${ECCKD_PREFIX}_lw_${INCODE}_${MODEL_CODE}.nc
	OUTPUT=${OUTDIR}/${ECCKD_PREFIX}_lw_${OUTCODE}_${MODEL_CODE}.nc
	LOG=${OUTDIR}/${ECCKD_PREFIX}_lw_${OUTCODE}_${MODEL_CODE}.log

	$OPTIMIZE_LUT \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    gpointfile=${GPOINTFILE} \
	    model_id=lw_${APPLICATION}_${BANDSTRUCT}-tol${TOL} \
	    $EXTRA_ARGS \
	    config_optimize_lut_${OPTIMIZE_MODE}.cfg \
	    |& tee $LOG
	test "${PIPESTATUS[0]}" -eq 0

    done
done

done
