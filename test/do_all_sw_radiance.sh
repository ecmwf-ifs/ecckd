#!/bin/bash
# Master script for creating shortwave CKD models

. config.h

# 0. Settings
APPLICATION=climate
BAND_STRUCTURE="modis"
TOLERANCE="0.5"
APP=climate-linear

# Nominal central wavelength (um): 0.412, 0.47, 0.55, 0.65, 0.86, 1.24, 1.63 and 2.11 
# Wavelength bound 1 (nm): 405 457 545 620 841 1230 1628 2105 
# Wavelength bound 2 (nm): 420 479 565 670 876 1250 1652 2155 
# Wavenumber bounds (cm-1):
WN1_SW_CUSTOM="23810 20877 17699 14925 11416 8000 6053 4640"
WN2_SW_CUSTOM="24691 21882 18349 16129 11891 8130 6143 4751"

# Make variables available to scripts find_g_points_sw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE
export WN1_SW_CUSTOM
export WN2_SW_CUSTOM
export APP

# 1. Merge well-mixed gases
#./merge_well_mixed_sw.sh

# 2. Reorder spectra
#./reorder_spectrum_sw.sh

# 3. Find g-points
#./find_g_points_sw.sh

# 4. Create raw CKD look-up table
./create_lut_sw.sh

# 4b. Scale look-up table to be exact for median profile
#./scale_lut_sw.sh

# 5. Optimize CKD look-up table
#./optimize_lut_sw.sh $OPTIMIZE_MODE_LIST

# 6. Run two-stream radiative transfer or just compute optical depths
# for CKDMIP scenarios
#./run_ckd_sw.sh
