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
# Master script for creating shortwave CKD models

. config.h

# 0. Settings

# First select the "application"
#APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
APPLICATION=climate

# Optionally specify a suffix to use for a particular experiment
MODEL_CODE_SUFFIX=
#MODEL_CODE_SUFFIX=-caviar

# Optionally select alternative water vapour continuum model
H2OCONTINUUM=
#H2OCONTINUUM=caviar

if [ "$APPLICATION" = climate ]
then
    # The best strategy for climate CKD models is to optimize first
    # H2O, CO2, O3 and O2+N2, then in subsequent optimization steps to
    # do CH4 and N2O. The CFCs may be ignored in the shortwave.
    OPTIMIZE_MODE_LIST="relative-base relative-ch4 relative-n2o"
else
    unset OPTIMIZE_MODE_LIST
fi

# Create 14 CKD models
BAND_STRUCTURE="wide narrow"
TOLERANCE="0.6 0.4 0.2 0.15 0.1 0.05 0.025"

# Create just one CKD models
BAND_STRUCTURE="wide"
TOLERANCE="0.2"

# "Red-green-blue" band structure with 16 and/or 32 g-points
BAND_STRUCTURE="rgb"
#TOLERANCE="0.16 0.047"
TOLERANCE="0.047"

# Or with these number of g-points: 12 16 20 24 28 32 36 40 48 64.
#TOLERANCE="0.3 0.16 0.11 0.072 0.062 0.047 0.04 0.03 0.0235 0.0121"

# "fine" band structure for reference calculations
#BAND_STRUCTURE="fine"
#TOLERANCE=0.0302

#BAND_STRUCTURE="double"
#TOLERANCE="0.065"

# "Reference" model with 64 points
#BAND_STRUCTURE=narrow
#TOLERANCE=0.019

# Another 64-point model with bands organised around the near-IR
# windows plus fine structure in the UV for UV index calculations
#BAND_STRUCTURE=window
#TOLERANCE=0.020

# Very fine band structure for diagnostics (96 g-points)
#BAND_STRUCTURE="vfine"
#TOLERANCE=0.02

# H2O suffix to accommodate different continuum models
if [ ! "$H2OCONTINUUM" ]
then
    # Default continuum
    H2OSUFFIX=
else
    H2OSUFFIX=-$H2OCONTINUUM
fi

# Make variables available to scripts find_g_points_sw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE
export MODEL_CODE_SUFFIX
export H2OSUFFIX

# 1. Merge well-mixed gases
./merge_well_mixed_sw.sh

# 2. Reorder spectra
./reorder_spectrum_sw.sh

# 3. Find g-points
./find_g_points_sw.sh

# 4. Create raw CKD look-up table
./create_lut_sw.sh

# 4b. Scale look-up table to be exact for median profile
./scale_lut_sw.sh

# 5. Optimize CKD look-up table
./optimize_lut_sw.sh $OPTIMIZE_MODE_LIST

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios; first optionally specify which CKD models to
# run radiative transfer on (default is the final one, i.e. "ckd", but
# "scaled" and "raw" are also possible)
#export VERSIONS="ckd"
./run_ckd_sw.sh

# 7. Copy to the CKDMIP directory, changing file names from stating
# tolerance to stating the total number of g points
#./copy_to_ckdmip_sw.sh
