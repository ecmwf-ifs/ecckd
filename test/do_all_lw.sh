#!/bin/bash
# Master script for creating longwave CKD models

. config.h

# 0. Settings
#APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
APPLICATION=climate

BAND_STRUCTURE="fsck wide narrow"
TOLERANCE="0.16 0.08 0.04 0.02 0.01 0.005"
TOLERANCE="0.04 0.01"

# Make variables available to scripts find_g_points_lw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE

# 1. Merge well-mixed gases
./merge_well_mixed_lw.sh

# 2. Reorder spectra
./reorder_spectrum_lw.sh

# 3. Find g-points
#./find_g_points_lw.sh

# 4. Create raw CKD look-up table
./create_lut_lw.sh

# 5. Optimize CKD look-up table
./optimize_lut_lw.sh

#APPLICATION=climate4 MODEL_CODE_SUFFIX=-sep-rel ./optimize_lut_lw.sh
#APPLICATION=climate2 MODEL_CODE_SUFFIX=-sep-rel ./optimize_lut_lw.sh
#APPLICATION=climate3 MODEL_CODE_SUFFIX=-sep ./optimize_lut_lw.sh
#APPLICATION=climate MODEL_CODE_SUFFIX=-sep ./optimize_lut_lw.sh
#APPLICATION=climate2 MODEL_CODE_SUFFIX=-sep ./optimize_lut_lw.sh

# 5.1 Second step in the climate case
if [ "$APPLICATION" = climate ]
then
    APPLICATION=climate2 ./optimize_lut_lw.sh
fi

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios
#MODEL_CODE_SUFFIX=-sep-rel ./run_lw_ckd.sh

./run_lw_ckd.sh
