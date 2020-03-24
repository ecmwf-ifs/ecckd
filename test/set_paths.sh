# This script is included by the others in this directory, and sets
# the locations of data and executables on your system

# CKDMIP installation directory and CKDMIP executables
CKDMIP_DIR=/home/pa/parr/src/ckdmip-0.9
CKDMIP_TOOL=${CKDMIP_DIR}/bin/ckdmip_tool
CKDMIP_LW=${CKDMIP_DIR}/bin/ckdmip_lw
CKDMIP_SW=${CKDMIP_DIR}/bin/ckdmip_sw

# Directory for ecCKD executables
BINDIR=../src/ecckd
REORDER_SPECTRUM_LW=${BINDIR}/reorder_spectrum_lw
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

# Idealized ataset
IDEALIZED_CODE=idealized
IDEALIZED_DIR=${CKDMIP_DATA_DIR}/${IDEAL_CODE}
IDEALIZED_CONC_DIR=${IDEALIZED_DIR}/conc
IDEALIZED_LW_SPECTRA_DIR=${IDEALIZED_DIR}/lw_spectra

# Training and evaluation dataset
TRAINING_CODE=evaluation1
TRAINING_DIR=${CKDMIP_DATA_DIR}/${TRAINING_CODE}
TRAINING_CONC_DIR=${TRAINING_DIR}/conc
TRAINING_LW_SPECTRA_DIR=${TRAINING_DIR}/lw_spectra

# Work directory
#WORK_DIR=${SCRATCH}/fsck
WORK_DIR=/hugetmp/parr/fsck
WORK_LW_SPECTRA_DIR=${WORK_DIR}/lw_spectra
WORK_LW_ORDER_DIR=${WORK_DIR}/lw_order
WORK_LW_GPOINTS_DIR=${WORK_DIR}/lw_gpoints
WELL_MIXED_LW_SPECTRA=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_merge-well-mixed_present.h5
WELL_MIXED_LW_SPECTRA_MINIMUM=${WORK_LW_SPECTRA_DIR}/ckdmip_${MMM_CODE}_lw_spectra_merge-well-mixed_minimum.h5

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
