#!/bin/bash
# Master script for creating longwave CKD models

. config.h

# 0. Settings
APPLICATION=global-nwp
BAND_STRUCTURE=microwave
APP=nwp-microwave
TOLERANCE=0.5

# Nominal central frequency (GHz): 31 166
# Bandpass (GHz): 0.2 2
WN1_LW_CUSTOM="1.03071 5.50381"
WN2_LW_CUSTOM="1.03738 5.57052"

# Make variables available to scripts find_g_points_lw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE
export WN1_LW_CUSTOM
export WN2_LW_CUSTOM
export APP

# 1. Merge well-mixed gases
#./merge_well_mixed_lw.sh

# 2. Reorder spectra
#./reorder_spectrum_lw.sh

# 3. Find g-points
./find_g_points_lw.sh

# 4. Create raw CKD look-up table
./create_lut_lw.sh

# 5. Optimize CKD look-up table
#./optimize_lut_lw.sh $OPTIMIZE_MODE_LIST

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios
#./run_ckd_lw.sh
