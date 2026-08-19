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
#APP=nwp-microwave
TOLERANCE=0.5

#APPLICATION=global-nwp
#BAND_STRUCTURE=msi

# Nominal central frequency (GHz): 31 166
# Bandpass (GHz): 0.2 2.8
#WN1_LW_CUSTOM="1.03071 5.47379"
#WN2_LW_CUSTOM="1.03738 5.60054"

# MSI thermal channels
#WN1_LW_CUSTOM="1084 885 800"
#WN2_LW_CUSTOM="1195 976 870"

# MODIS thermal channels
#6.535–6.895
#8.400–8.700
#10.780–11.280
#11.770–12.270

#BAND_STRUCTURE=modis
#WN1_LW_CUSTOM="1450 1149 887 815"
#WN2_LW_CUSTOM="1530 1190 928 850"

# GMI microwave channels
#Channel No	Central Frequency (Ghz)	Central Frequency Stabilization (±MHz)	Bandwidth (Mhz)	Polarization	Integration time (ms)	NEDT (K)	Antenna beamwidth @ 3dB (º)
#1	10.65	10	100	V	9.7	0.96	1.75
#2	10.65	10	100	H	9.7	0.96	1.75
#3	18.70	20	200	V	5.3	0.84	1.00
#4	18.70	20	200	H	5.3	0.84	1.00
#5	23.80	20	400	V	5.0	1.05	0.90
#6	36.50	50	1000	V	5.0	0.65	0.90
#7	36.5	50	1000	H	5.0	0.65	0.90
#8	89.00	200	6000	V	2.2	0.57	0.40
#9	89.00	200	6000	H	2.2	0.57	0.40
#10	166.0	200	3000	V	3.6	1.5	0.40
#11	166.0	200	3000	H	3.6	1.5	0.40
#12	183.31±3	200	3500	V	3.6	1.5	0.4
#13	183.31±7	200	4500	V	3.6	1.5	0.4

BAND_STRUCTURE=gmi
WN1_LW_CUSTOM=""
WN2_LW_CUSTOM=""

# Make variables available to scripts find_g_points_lw.sh onwards
export TOLERANCE
export APPLICATION
export BAND_STRUCTURE
export WN1_LW_CUSTOM
export WN2_LW_CUSTOM
#export APP

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
# for CKDMIP scenarios
#./run_ckd_lw.sh
