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
# Master script for creating longwave gas-optics models for microwave
# remote sensing; note that these models are monochromatic so no
# optimization step is needed.

. config.h

# 0. Settings
APPLICATION=global-nwp
BAND_STRUCTURE=microwave
APP=nwp-microwave
TOLERANCE=0.5

# Nominal central frequency (GHz): 31 166
# Bandpass (GHz): 0.2 2.8
WN1_LW_CUSTOM="1.03071 5.47379"
WN2_LW_CUSTOM="1.03738 5.60054"

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
