#!/bin/bash
# Run two-stream radiative transfer. The input requirements are the same as
# find_g_points_sw.sh 

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

# If the version of the model to run is not stated, do just the final
# one
if [ -z "$VERSIONS" ]
then
    VERSIONS=ckd
fi

# The only RT option in the shortwave is two-stream, denoted by
# "fluxes"
FLUXESSTR=fluxes

# Which scenarios need to be simulated
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

cat > config_sw_ckd_rt_evaluation.nam <<EOF
&shortwave_config
optical_depth_name       = "optical_depth",
rayleigh_optical_depth_name = "rayleigh_optical_depth",
pressure_name            = "pressure_hl",
incoming_flux_name = "incoming_sw",
do_write_spectral_fluxes = true,
do_write_optical_depth   = false,
surf_albedo = 0.15,
use_mu0_dimension = true,
cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
iverbose = 2
/
EOF

mkdir -p ${WORK_SW_CKD_OD_DIR}
mkdir -p ${WORK_SW_FLUX_DIR}

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
	    CKD_MODEL=${WORK_DIR}/sw_${VERDIR}/${ECCKD_PREFIX}_sw_${VERSTR}_${MODEL_CODE}.nc

	    ${BANNER} Running CKD model ${MODEL_CODE} for relevant scenarios

	    for SCENARIO in $SCENARIOS
	    do
		INPUT=${EVALUATION_CONC_DIR}/ckdmip_${EVALUATION_CODE}_concentrations_${SCENARIO}.nc
		EVALUATION_PREFIX=${ECCKD_PREFIX}_${EVALUATION_CODE}_sw_${MODEL_CODE}${VERSUFFIX}
		OD_FILE=${WORK_SW_CKD_OD_DIR}/${EVALUATION_PREFIX}_optical-depth_${SCENARIO}.nc
		FLUX_FILE=$WORK_SW_FLUX_DIR/${EVALUATION_PREFIX}_${FLUXESSTR}_${SCENARIO}.nc

		# Compute spectral optical depths from CKD model
		${BANNER_MINI} Computing optical depths for scenario $SCENARIO
		$SW_CKD ckd_model=$CKD_MODEL input=${INPUT} output=${OD_FILE} tsi=1361.0
# write_od_only=1

		# Compute fluxes from optical depths
		${BANNER_MINI} Running radiative transfer for scenario $SCENARIO
		$CKDMIP_SW --config config_sw_ckd_rt_evaluation.nam \
		    --scenario $SCENARIO \
		    --ckd $OD_FILE \
		    --output $FLUX_FILE
	    done
	done
    done
done
