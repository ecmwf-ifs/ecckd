# -*- shell-script -*-
# This script is included by the others in this directory, and sets
# the locations of data and executables on your system

# Fail on error
set -e

# Use all threads available
unset OMP_NUM_THREADS

# CKDMIP installation directory and CKDMIP executables
CKDMIP_DIR=/home/pa/parr/src/ckdmip-0.9
CKDMIP_TOOL=${CKDMIP_DIR}/bin/ckdmip_tool
CKDMIP_LW=${CKDMIP_DIR}/bin/ckdmip_lw
CKDMIP_SW=${CKDMIP_DIR}/bin/ckdmip_sw

# Directory for ecCKD executables
BINDIR=../src/ecckd
REORDER_SPECTRUM="${BINDIR}/reorder_spectrum"
FIND_G_POINTS_LW=${BINDIR}/find_g_points_lw
CREATE_LOOK_UP_TABLE=${BINDIR}/create_look_up_table
OPTIMIZE_LUT=${BINDIR}/optimize_lut
LW_SPECTRA=${BINDIR}/lw_spectra
LW_CKD=${BINDIR}/lw_ckd

# CKDMIP data directories
CKDMIP_DATA_DIR=/hugetmp/parr/ckdmip

# Median/minimum/maximum dataset
MMM_CODE=mmm
MMM_DIR=${CKDMIP_DATA_DIR}/${MMM_CODE}
MMM_CONC_DIR=${MMM_DIR}/conc
MMM_LW_SPECTRA_DIR=${MMM_DIR}/lw_spectra
MMM_SW_SPECTRA_DIR=${MMM_DIR}/sw_spectra
MMM_SW_SPECTRA_EXTRA_DIR=${MMM_DIR}/sw_spectra_extras
MMM_SW_SSI=${MMM_SW_SPECTRA_EXTRA_DIR}/ckdmip_ssi.h5

# Idealized ataset
IDEALIZED_CODE=idealized
IDEALIZED_DIR=${CKDMIP_DATA_DIR}/${IDEALIZED_CODE}
IDEALIZED_CONC_DIR=${IDEALIZED_DIR}/conc
IDEALIZED_LW_SPECTRA_DIR=${IDEALIZED_DIR}/lw_spectra
IDEALIZED_SW_SPECTRA_DIR=${IDEALIZED_DIR}/sw_spectra

# Training and evaluation dataset
TRAINING_CODE=evaluation1
TRAINING_DIR=${CKDMIP_DATA_DIR}/${TRAINING_CODE}
TRAINING_CONC_DIR=${TRAINING_DIR}/conc
TRAINING_LW_SPECTRA_DIR=${TRAINING_DIR}/lw_spectra
TRAINING_LW_FLUXES_DIR=${TRAINING_DIR}/lw_fluxes
TRAINING_SW_SPECTRA_DIR=${TRAINING_DIR}/sw_spectra
TRAINING_SW_FLUXES_DIR=${TRAINING_DIR}/sw_fluxes

# Work directory
WORK_DIR=/hugetmp/parr/ecckd
WORK_LW_SPECTRA_DIR=${WORK_DIR}/lw_spectra
WORK_LW_ORDER_DIR=${WORK_DIR}/lw_order
WORK_LW_GPOINTS_DIR=${WORK_DIR}/lw_gpoints
WORK_LW_RAW_CKD_DIR=${WORK_DIR}/lw_raw-ckd
WORK_LW_CKD_DIR=${WORK_DIR}/lw_ckd
WORK_LW_CKD_OD_DIR=${WORK_DIR}/lw_optical-depth
WORK_SW_SPECTRA_DIR=${WORK_DIR}/sw_spectra
WORK_SW_ORDER_DIR=${WORK_DIR}/sw_order
WORK_SW_GPOINTS_DIR=${WORK_DIR}/sw_gpoints
WORK_SW_RAW_CKD_DIR=${WORK_DIR}/sw_raw-ckd
WORK_SW_CKD_DIR=${WORK_DIR}/sw_ckd
WORK_SW_CKD_OD_DIR=${WORK_DIR}/sw_optical-depth

WELL_MIXED_LW_SPECTRA=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_composite_present.h5
WELL_MIXED_LW_SPECTRA_MINIMUM=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_composite_minimum.h5
WELL_MIXED_SW_SPECTRA=${WORK_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_composite_present.h5
WELL_MIXED_SW_SPECTRA_MINIMUM=${WORK_SW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_sw_spectra_composite_minimum.h5

# Band definitions from CKDMIP
WN1_LW_NARROW="0 350 500 630 700 820 980 1080 1180 1390 1480 1800 2080"
WN2_LW_NARROW="350 500 630 700 820 980 1080 1180 1390 1480 1800 2080 3260"
WN1_LW_WIDE="0 500 820 1180 1800"
WN2_LW_WIDE="500 820 1180 1800 3260"
WN1_SW_NARROW="250 2600 3250 4000 4650 5150 6150 8050 12850 16000 22650 29000 38000"
WN2_SW_NARROW="2600 3250 4000 4650 5150 6150 8050 12850 16000 22650 29000 38000 50000"
WN1_SW_WIDE="250 4000 8050 16000 29000"
WN2_SW_WIDE="4000 8050 16000 29000 50000"

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
