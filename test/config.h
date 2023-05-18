# -*- shell-script -*-
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
# This script is included by the others in this directory, and sets
# the locations of data and executables on your system

# Fail on error
set -e

# Use all threads available
unset OMP_NUM_THREADS

# Include file defining ecCKD version
. version.h

# CKDMIP installation directory and CKDMIP executables
CKDMIP_DIR=/home/parr/git/ckdmip
CKDMIP_TOOL=${CKDMIP_DIR}/bin/ckdmip_tool
CKDMIP_LW=${CKDMIP_DIR}/bin/ckdmip_lw
CKDMIP_SW=${CKDMIP_DIR}/bin/ckdmip_sw

# Directory for ecCKD executables
BINDIR=../src/ecckd
REORDER_SPECTRUM="${BINDIR}/reorder_spectrum"
REORDER_CLOUD_SPECTRUM="${BINDIR}/reorder_cloud_spectrum"
FIND_G_POINTS=${BINDIR}/find_g_points
CREATE_LOOK_UP_TABLE=${BINDIR}/create_look_up_table
OPTIMIZE_LUT=${BINDIR}/optimize_lut
SCALE_LUT=${BINDIR}/scale_lut
LW_SPECTRA=${BINDIR}/lw_spectra
LW_CKD=${BINDIR}/run_ckd
SW_CKD=${BINDIR}/run_ckd

# CKDMIP data directories
CKDMIP_DATA_DIR=/perm/parr/ckdmip

# Median/minimum/maximum dataset
MMM_CODE=mmm
MMM_DIR=${CKDMIP_DATA_DIR}/${MMM_CODE}
MMM_CONC_DIR=${MMM_DIR}/conc
MMM_LW_SPECTRA_DIR=${MMM_DIR}/lw_spectra
MMM_SW_SPECTRA_DIR=${MMM_DIR}/sw_spectra
MMM_SW_SPECTRA_EXTRA_DIR=${MMM_DIR}/sw_spectra_extras
MMM_SW_SSI=${MMM_SW_SPECTRA_EXTRA_DIR}/ckdmip_ssi.h5

# Idealized dataset
IDEALIZED_CODE=idealized
IDEALIZED_DIR=${CKDMIP_DATA_DIR}/${IDEALIZED_CODE}
IDEALIZED_CONC_DIR=${IDEALIZED_DIR}/conc
IDEALIZED_LW_SPECTRA_DIR=${IDEALIZED_DIR}/lw_spectra
IDEALIZED_SW_SPECTRA_DIR=${IDEALIZED_DIR}/sw_spectra

# Training and evaluation dataset (usually evaluation1)
TRAINING_CODE=evaluation1
#TRAINING_CODE=evaluation1-caviar
#TRAINING_CODE=evaluation1-cropsurf2
TRAINING_DIR=${CKDMIP_DATA_DIR}/${TRAINING_CODE}
TRAINING_CONC_DIR=${TRAINING_DIR}/conc
TRAINING_LW_SPECTRA_DIR=${TRAINING_DIR}/lw_spectra
TRAINING_LW_FLUXES_DIR=${TRAINING_DIR}/lw_fluxes
TRAINING_SW_SPECTRA_DIR=${TRAINING_DIR}/sw_spectra
TRAINING_SW_FLUXES_DIR=${TRAINING_DIR}/sw_fluxes

#TRAINING_SW_SSI=${TRAINING_SW_SPECTRA_DIR}/ckdmip_ssi.h5
TRAINING_SW_SSI=${CKDMIP_DATA_DIR}/evaluation1/sw_spectra/ckdmip_ssi.h5

# Do we train with both evaluation1 and evaluation2 ("yes" or "no")
TRAINING_BOTH=no
# When we train with "both" what is the code of the second training
# dataset?
TRAINING_CODE2=evaluation2

# Evaluation dataset (could be evaluation1 or evaluation2)
EVALUATION_CODE=evaluation2
EVALUATION_DIR=${CKDMIP_DATA_DIR}/${EVALUATION_CODE}
EVALUATION_CONC_DIR=${EVALUATION_DIR}/conc
EVALUATION_LW_SPECTRA_DIR=${EVALUATION_DIR}/lw_spectra
EVALUATION_LW_FLUXES_DIR=${EVALUATION_DIR}/lw_fluxes
EVALUATION_SW_SPECTRA_DIR=${EVALUATION_DIR}/sw_spectra
EVALUATION_SW_FLUXES_DIR=${EVALUATION_DIR}/sw_fluxes

# Cloud spectrum
CLOUD_SPECTRUM=../data/mie_droplet_scattering.nc

# Work directory
WORK_DIR=/perm/parr/ecckd
WORK_LW_SPECTRA_DIR=${WORK_DIR}/lw_spectra
WORK_LW_ORDER_DIR=${WORK_DIR}/lw_order
WORK_LW_GPOINTS_DIR=${WORK_DIR}/lw_gpoints
WORK_LW_RAW_CKD_DIR=${WORK_DIR}/lw_raw-ckd-definition
WORK_LW_CKD_DIR=${WORK_DIR}/lw_ckd-definition
WORK_LW_CKD_OD_DIR=${WORK_DIR}/lw_optical-depth
WORK_LW_LBL_FLUX_DIR=${WORK_DIR}/lw_lbl_fluxes
WORK_LW_FLUX_DIR=${WORK_DIR}/lw_fluxes
WORK_SW_SPECTRA_DIR=${WORK_DIR}/sw_spectra
WORK_SW_ORDER_DIR=${WORK_DIR}/sw_order
WORK_SW_GPOINTS_DIR=${WORK_DIR}/sw_gpoints
WORK_SW_RAW_CKD_DIR=${WORK_DIR}/sw_raw-ckd-definition
WORK_SW_CKD_DIR=${WORK_DIR}/sw_ckd-definition
WORK_SW_CKD_OD_DIR=${WORK_DIR}/sw_optical-depth
WORK_SW_LBL_FLUX_DIR=${WORK_DIR}/sw_lbl_fluxes
WORK_SW_FLUX_DIR=${WORK_DIR}/sw_fluxes

WELL_MIXED_LW_SPECTRA=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_composite_present.h5
WELL_MIXED_LW_SPECTRA_MINIMUM=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_composite_minimum.h5
WELL_MIXED_LW_SPECTRA_O2N2=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_o2n2_constant.h5

WELL_MIXED_SW_SPECTRA=${WORK_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_composite_present.h5
WELL_MIXED_SW_SPECTRA_MINIMUM=${WORK_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_composite_minimum.h5
WELL_MIXED_SW_SPECTRA_O2N2=${WORK_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_o2n2_constant.h5

# Band definitions from CKDMIP
WN1_LW_NARROW="0 350 500 630 700 820 980 1080 1180 1390 1480 1800 2080"
WN2_LW_NARROW="350 500 630 700 820 980 1080 1180 1390 1480 1800 2080 3260"
WN1_LW_WIDE="0 500 820 1180 1800"
WN2_LW_WIDE="500 820 1180 1800 3260"
WN1_SW_NARROW="250 2600 3250 4000 4650 5150 6150 8050 12850 16000 22650 29000 38000"
WN2_SW_NARROW="2600 3250 4000 4650 5150 6150 8050 12850 16000 22650 29000 38000 50000"
WN1_SW_WIDE="250 4000 8050 16000 29000"
WN2_SW_WIDE="4000 8050 16000 29000 50000"
WN1_SW_DOUBLE="250 16000"
WN2_SW_DOUBLE="16000 50000"
WN1_SW_RGB="250 14300 16650 20000 25000"
WN2_SW_RGB="14300 16650 20000 25000 50000"
WN1_SW_GB="250 8000 16650 20000 25000"
WN2_SW_GB="8000 16650 20000 25000 50000"
WN1_SW_FINE="250 3750 5350 7150 8700 10650 12100 13350 14300 15400 16650 18200 20000 22200 25000 28550 30250 30750 31250 31750 32250 32750 33250 33750 34250"
WN2_SW_FINE="3750 5350 7150 8700 10650 12100 13350 14300 15400 16650 18200 20000 22200 25000 28550 30250 30750 31250 31750 32250 32750 33250 33750 34250 50000"
WN1_SW_WINDOW="250 3750 5350 7150 8700 10650 14300 16650 20000 25000 28550 30250 30750 31250 31750 32250 32750 33250 33750"
WN2_SW_WINDOW="3750 5350 7150 8700 10650 14300 16650 20000 25000 28550 30250 30750 31250 31750 32250 32750 33250 33750 50000"
WN1_SW_VFINE="250 2600 3750 5350 7150 8700 10650 12100 13350 13800 14300 14800 15400 16000 16650 17400 18200 19050 20000 21050 22200 23550 25000 26300 26650 27050 27400 27800 28150 28550 29000 29400 29850 30300 30750 31250 31750 32250 32800 33350 33900 34500 35100 35700"
WN2_SW_VFINE="2600 3750 5350 7150 8700 10650 12100 13350 13800 14300 14800 15400 16000 16650 17400 18200 19050 20000 21050 22200 23550 25000 26300 26650 27050 27400 27800 28150 28550 29000 29400 29850 30300 30750 31250 31750 32250 32800 33350 33900 34500 35100 35700 50000"

# Prefix final files by the following
ECCKD_PREFIX=ecckd-$ECCKD_VERSION

CKDMIP_RESULTS_DIR=${CKDMIP_DATA_DIR}/results/${ECCKD_PREFIX}

# Should the shortwave composites include Rayleigh scattering for
# ordering the spectra?
COMPOSITE_SW_INCLUDES_RAYLEIGH=yes

function my_banner {
    echo
    echo "####################################################################"
    echo "### ${@^^}"
    echo "####################################################################"
}

function my_banner_skip {
    echo
    echo "### ${@}"
}

function my_banner_error {
    echo "********************************************************************"
    echo "*** ${@^^}"
    echo "********************************************************************"
}

BANNER=my_banner
BANNER_SKIP=my_banner_skip
BANNER_ERROR=my_banner_skip
BANNER_MINI=my_banner_skip
