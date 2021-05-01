#!/bin/bash
# Master script for creating shortwave CKD models

. config.h

# 0. Settings
APPLICATION=limited-area-nwp
APPLICATION=global-nwp
APPLICATION=climate

if [ "$APPLICATION" = climate ]
then
    # The best strategy for climate CKD models is to optimize first
    # H2O, CO2, O3 and O2+N2, then in subsequent optimization steps to
    # do CH4 and N2O. The CFCs may be ignored in the shortwave.
    OPTIMIZE_MODE_LIST="relative-base relative-ch4 relative-n2o"
    #OPTIMIZE_MODE_LIST="relative-base relative-minor" # Less good than CH4 and N2O separate
    #OPTIMIZE_MODE_LIST="relative-ch4 relative-n2o"
    #OPTIMIZE_MODE_LIST="relative-base"
else
    unset OPTIMIZE_MODE_LIST
fi

#BAND_STRUCTURE="fsck double wide narrow"
BAND_STRUCTURE="wide narrow"
TOLERANCE="0.6 0.4 0.2 0.15 0.1 0.05 0.025"

BAND_STRUCTURE="wide"
TOLERANCE="0.2 0.1 0.05 0.025"

BAND_STRUCTURE="rgb"
TOLERANCE="0.4 0.2 0.8"
TOLERANCE="3.0 5.0"
# Make variables available to scripts find_g_points_sw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE

# 1. Merge well-mixed gases
#./merge_well_mixed_sw.sh

# 2. Reorder spectra
./reorder_spectrum_sw.sh

# 3. Find g-points
#./find_g_points_sw.sh

# 4. Create raw CKD look-up table
#./create_lut_sw.sh

# 4b. Scale look-up table to be exact for median profile
#./scale_lut_sw.sh

# 5. Optimize CKD look-up table
#./optimize_lut_sw.sh $OPTIMIZE_MODE_LIST

export VERSIONS="ckd"

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios
#./run_ckd_sw.sh

# 7. Copy
#./copy_to_ckdmip_sw.sh
