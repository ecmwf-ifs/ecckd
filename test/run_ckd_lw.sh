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
# This script runs two-stream longwave radiative transfer. The input
# requirements are the same as find_g_points_lw.sh. It is called from
# do_all_lw.sh and should not be run directly.

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

# If the version of the model to run is not stated, do just the final
# one
if [ -z "$VERSIONS" ]
then
    VERSIONS=ckd
fi

# Number of angles per hemisphere
NANGLE=4

if [ "$NANGLE" = 0 ]
then
    FLUXESSTR=fluxes
else
    FLUXESSTR=fluxes-${NANGLE}angle
fi

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
		$LW_CKD ckd_model=$CKD_MODEL input=${INPUT} output=${OD_FILE}
# write_od_only=1

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
