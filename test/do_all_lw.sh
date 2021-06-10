#!/bin/bash
# Master script for creating longwave CKD models

. config.h

# 0. Settings

# First select the "application"
#APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
APPLICATION=climate

if [ "$APPLICATION" = climate ]
then
    # The best strategy for climate CKD models is to optimize first
    # H2O, CO2, O3 and O2+N2, then in subsequent optimization steps to
    # do the three next minor gases.  CFC12 is already good enough
    # without optimization.
    OPTIMIZE_MODE_LIST="relative-base relative-ch4 relative-n2o relative-cfc11"
else
    unset OPTIMIZE_MODE_LIST
fi

# Create 18 CKD models
BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.16 0.08 0.04 0.02 0.01 0.005"

# Create just one CKD model
BAND_STRUCTURE=fsck
TOLERANCE=0.02

# Make variables available to scripts find_g_points_lw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE

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
#export VERSIONS="ckd"
./run_ckd_lw.sh
