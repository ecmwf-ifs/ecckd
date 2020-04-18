#!/bin/bash
# Master script for creating CKD models

. config.h

# 0. Settings
#APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
APPLICATION=climate

BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.16 0.08 0.04 0.02 0.01 0.005"

# Make variables available to scripts find_g_points_lw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE

# 1. Merge well-mixed gases
./merge_well_mixed.sh

# 2. Reorder spectra
./reorder_spectrum_lw.sh

# 3. Find g-points
./find_g_points_lw.sh

# 4. Create raw CKD look-up table
./create_lut_lw.sh

# 5. Optimize CKD look-up table
./optimize_lut_lw.sh

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios
./run_lw_ckd.sh
