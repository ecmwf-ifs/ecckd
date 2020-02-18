#!/bin/bash
# Master script for creating CKD models

set -e

. set_paths.sh

# 0. Settings
#APPLICATION=limited-area-nwp
#APPLICATION=global-nwp
APPLICATION=climate

CONFIGURATION=fsck
TOLERANCE=0.01

# Consolidate settings

# Select minimum pressure to quantify heating-rate errors (Pa)
if [ ${APPLICATION} = limited-area-nwp ]
then
    MIN_PRESSURE=400
    APP=nwp
else
    MIN_PRESSURE=2
    if [ ${APPLICATION} = climate ]
    then
	APP=climate
    else
	APP=nwp
    fi
fi

MODEL_CODE=${APPLICATION}_${CONFIGURATION}-tol${TOLERANCE}

# 1. Merge well-mixed gases

./merge_well_mixed.sh

# 2. Reorder spectra

./reorder_spectrum_lw.sh

# 3. Find g-points

MODEL_CODE=${MODEL_CODE} MIN_PRESSURE=${MIN_PRESSURE} TOLERANCE=${TOLERANCE} APP=${APP} ./find_g_points_lw.sh

