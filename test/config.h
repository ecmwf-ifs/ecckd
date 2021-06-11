# -*- shell-script -*-
# This script is included by the others in this directory, and sets
# the locations of data and executables on your system

# Fail on error
set -e

# Use all threads available
unset OMP_NUM_THREADS

# Include file defining ecCKD version
. version.h

# CKDMIP installation directory and CKDMIP executables
CKDMIP_DIR=/home/pa/parr/src/ckdmip-1.0
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
CKDMIP_DATA_DIR=/hugetmp/parr/ckdmip

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
TRAINING_DIR=${CKDMIP_DATA_DIR}/${TRAINING_CODE}
TRAINING_CONC_DIR=${TRAINING_DIR}/conc
TRAINING_LW_SPECTRA_DIR=${TRAINING_DIR}/lw_spectra
TRAINING_LW_FLUXES_DIR=${TRAINING_DIR}/lw_fluxes
TRAINING_SW_SPECTRA_DIR=${TRAINING_DIR}/sw_spectra
TRAINING_SW_FLUXES_DIR=${TRAINING_DIR}/sw_fluxes
TRAINING_SW_SSI=${TRAINING_SW_SPECTRA_DIR}/ckdmip_ssi.h5

# Evaluation dataset (could be evaluation1 or evaluation2)
EVALUATION_CODE=evaluation1
EVALUATION_DIR=${CKDMIP_DATA_DIR}/${EVALUATION_CODE}
EVALUATION_CONC_DIR=${EVALUATION_DIR}/conc
EVALUATION_LW_SPECTRA_DIR=${EVALUATION_DIR}/lw_spectra
EVALUATION_LW_FLUXES_DIR=${EVALUATION_DIR}/lw_fluxes
EVALUATION_SW_SPECTRA_DIR=${EVALUATION_DIR}/sw_spectra
EVALUATION_SW_FLUXES_DIR=${EVALUATION_DIR}/sw_fluxes

# Cloud spectrum
CLOUD_SPECTRUM=../data/mie_droplet_scattering.nc

# Work directory
WORK_DIR=/hugetmp/parr/ecckd
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
WN1_SW_DOUBLE="250 12850"
WN2_SW_DOUBLE="12850 50000"
WN1_SW_RGB="250 14300 16650 20000 25000"
WN2_SW_RGB="14300 16650 20000 25000 50000"
WN1_SW_GB="250 8000 16650 20000 25000"
WN2_SW_GB="8000 16650 20000 25000 50000"
#WN1_SW_GB="250 16650 20000 25000"
#WN2_SW_GB="16650 20000 25000 50000"

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
