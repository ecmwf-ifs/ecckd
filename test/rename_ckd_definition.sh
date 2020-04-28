#!/bin/bash
. config.h

set -ex

VERSIONS="ckd"
APPLICATION=climate
BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.04 0.01"
#MODEL_CODE_SUFFIX=-sep

mkdir -p ${WORK_DIR}/lw_spectral-definition/

# Loop over each band structure, tolerance, version (raw or final) and
# scenario
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	for VER in $VERSIONS
	do
	    MODEL_CODE=${APPLICATION}_${BANDSTRUCT}_tol${TOL}${MODEL_CODE_SUFFIX}
	    cd ${WORK_DIR}/lw_spectral-definition/
	    ln -v -s ../lw_ckd-definition/lw_ckd-definition_${MODEL_CODE}.nc \
		ecckd_lw_${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}_spectral-definition.nc
	done
    done
done
