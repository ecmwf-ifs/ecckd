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
# Master script for creating longwave CKD models

. config.h

# 0. Settings

# First select the "application"
#APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
APPLICATION=climate

# Optionally specify a suffix to use for a particular experiment
#MODEL_CODE_SUFFIX=-no-continuum
MODEL_CODE_SUFFIX=

# Optionally select alternative water vapour continuum model
#H2OCONTINUUM=no-continuum
H2OCONTINUUM=

if [ "$APPLICATION" = climate ]
then
    # The best strategy for climate CKD models is to optimize first
    # H2O, CO2, O3 and O2+N2, then in subsequent optimization steps to
    # do the three next minor gases.
    OPTIMIZE_MODE_LIST="relative-base relative-ch4 relative-n2o relative-cfc"
else
    unset OPTIMIZE_MODE_LIST
fi

# Create 18 CKD models
#BAND_STRUCTURE="fsck wide narrow"
#TOLERANCE="0.16 0.08 0.04 0.02 0.01 0.005"

# FSCK models with 16 and/or 32 g-points
BAND_STRUCTURE=fsck
#TOLERANCE="0.061 0.0161"
TOLERANCE="0.0161"

# Or with these number of g-points: 12 16 20 24 28 32 36 40 48 64.
# Note that the 24-point model has a large heating-rate bias at 0.03
# hPa, which can be overcome by setting prior_error=4.0 rather than
# 8.0 in optimize_lut_lw.sh
#TOLERANCE="0.11 0.061 0.043 0.03 0.02 0.0161 0.013 0.0105 0.00732 0.0047"

# Create a reference CKD model with 64 points
#BAND_STRUCTURE=narrow
#TOLERANCE=0.013

# H2O suffix to accommodate different continuum models
if [ ! "$H2OCONTINUUM" ]
then
    # Default continuum
    H2OSUFFIX=
else
    H2OSUFFIX=-$H2OCONTINUUM
fi

# Make variables available to scripts find_g_points_lw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE
export MODEL_CODE_SUFFIX
export H2OSUFFIX

# 1. Merge well-mixed gases
./merge_well_mixed_lw.sh

# 2. Reorder spectra
./reorder_spectrum_lw.sh

# 3. Find g-points
./find_g_points_lw.sh

# 4. Create raw CKD look-up table
./create_lut_lw.sh

# 5. Optimize CKD look-up table
./optimize_lut_lw.sh $OPTIMIZE_MODE_LIST

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios; first optionally specify which CKD models to
# run radiative transfer on (default is the final one, i.e. "ckd", but
# "raw", "raw2" etc are also possible)
#export VERSIONS="ckd "
./run_ckd_lw.sh

# 7. Copy to the CKDMIP directory, changing file names from stating
# tolerance to stating the total number of g points
#./copy_to_ckdmip_lw.sh
