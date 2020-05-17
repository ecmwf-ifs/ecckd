#!/bin/bash
# Master script for creating shortwave CKD models

. config.h

# 0. Settings
APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
#APPLICATION=climate

#BAND_STRUCTURE="fsck double wide narrow"
BAND_STRUCTURE="wide narrow"
TOLERANCE="0.8 0.4 0.2 0.1 0.05 0.025"
TOLERANCE="0.4 0.2 0.1"

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

# 5. Optimize CKD look-up table
./optimize_lut_sw.sh

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios
./run_ckd_sw.sh
