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
# This script renames ckd definition files from using their tolerance
# in the filename to their number of g points

. config.h

VERSIONS="ckd"
APPLICATION=climate
BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.16 0.08 0.04 0.02 0.01"

#MODEL_CODE_SUFFIX=-sep

BAND_STRUCTURE=fsck
TOLERANCE=0.02

mkdir -p ${WORK_DIR}/lw_spectral-definition/

# Loop over each band structure, tolerance, version (raw or final) and
# scenario
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	for VER in $VERSIONS
	do
	    MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	    cd ${WORK_DIR}/lw_spectral-definition/
	    ln -v -s ../lw_ckd-definition/${ECCKD_PREFIX}_lw_ckd-definition_${MODEL_CODE}.nc \
		${ECCKD_PREFIX}_lw_${MODEL_CODE}_spectral-definition.nc
	done
    done
done
