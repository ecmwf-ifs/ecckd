#!/bin/bash
# Master script for creating shortwave CKD models

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
    # do CH4 and N2O. The CFCs may be ignored in the shortwave.
    OPTIMIZE_MODE_LIST="relative-base relative-ch4 relative-n2o"
    #OPTIMIZE_MODE_LIST="relative-base relative-minor" # Less good than CH4 and N2O separate
else
    unset OPTIMIZE_MODE_LIST
fi

# Create 14 CKD models
BAND_STRUCTURE="wide narrow"
TOLERANCE="0.6 0.4 0.2 0.15 0.1 0.05 0.025"

# Create just one CKD models
BAND_STRUCTURE="wide"
TOLERANCE="0.1"

# Experimental "red-green-blue" band structure
#BAND_STRUCTURE="rgb"
#TOLERANCE=" 1.5"

# Make variables available to scripts find_g_points_sw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE

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
