#!/bin/bash
# Run two-stream radiative transfer. The input requirements are the same as
# find_g_points_lw.sh 

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

#VERSIONS="raw raw2 raw3 raw4 ckd"
VERSIONS=ckd
NANGLE=4

if [ "$NANGLE" = 0 ]
then
    FLUXESSTR=fluxes
else
    FLUXESSTR=fluxes-${NANGLE}angle
fi

if [ "$APP" = nwp ]
then
    SCENARIOS=present
else
    SCENARIOS="glacialmax
preindustrial
present
future
co2-180
co2-280
co2-560
co2-1120
co2-2240
ch4-350
ch4-700
ch4-1200
ch4-2600
ch4-3500
n2o-190
n2o-270
n2o-405
n2o-540
cfc11-0
cfc11-2000
cfc12-0
cfc12-550
co2-180-ch4-350
co2-2240-ch4-350
co2-180-ch4-3500
co2-2240-ch4-3500
co2-180-n2o-190
co2-2240-n2o-190
co2-180-n2o-540
co2-2240-n2o-540
ch4-350-n2o-190
ch4-3500-n2o-190
ch4-350-n2o-540
ch4-3500-n2o-540"

fi

cat > config_lw_ckd_rt_evaluation.nam <<EOF
&longwave_config
optical_depth_name       = "optical_depth",
pressure_name            = "pressure_hl",
nangle                   = $NANGLE,
do_write_planck          = false,
do_write_spectral_fluxes = true,
do_write_optical_depth   = false,
input_planck_per_sterad  = false,
iverbose = 2
/
EOF

mkdir -p ${WORK_LW_CKD_OD_DIR}
mkdir -p ${WORK_LW_FLUX_DIR}

# Loop over each band structure, tolerance, version (raw or final) and
# scenario
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	for VER in $VERSIONS
	do
	    if [ "$VER" = ckd ]
	    then
		VERSTR="${VER}-definition"
		VERSUFFIX=
		VERDIR=$VERSTR
	    else
		VERSTR="${VER}-ckd-definition"
		VERSUFFIX=-${VER}
		VERDIR=raw-ckd-definition
	    fi

	    MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	    CKD_MODEL=${WORK_DIR}/lw_${VERDIR}/${ECCKD_PREFIX}_lw_${VERSTR}_${MODEL_CODE}.nc

	    ${BANNER} Running CKD model ${MODEL_CODE} for relevant scenarios

	    for SCENARIO in $SCENARIOS
	    do
		INPUT=${EVALUATION_CONC_DIR}/ckdmip_${EVALUATION_CODE}_concentrations_${SCENARIO}.nc
		EVALUATION_PREFIX=${ECCKD_PREFIX}_${EVALUATION_CODE}_lw_${MODEL_CODE}${VERSUFFIX}
		OD_FILE=${WORK_LW_CKD_OD_DIR}/${EVALUATION_PREFIX}_optical-depth_${SCENARIO}.nc
		FLUX_FILE=$WORK_LW_FLUX_DIR/${EVALUATION_PREFIX}_${FLUXESSTR}_${SCENARIO}.nc

		# Compute spectral optical depths from CKD model
		${BANNER_MINI} Computing optical depths for scenario $SCENARIO
		$LW_CKD ckd_model=$CKD_MODEL input=${INPUT} output=${OD_FILE} write_od_only=1

		# Compute fluxes from optical depths
		${BANNER_MINI} Running radiative transfer for scenario $SCENARIO
		$CKDMIP_LW --config config_lw_ckd_rt_evaluation.nam \
		    --scenario $SCENARIO \
		    --ckd $OD_FILE \
		    --output $FLUX_FILE
	    done
	done
    done
done