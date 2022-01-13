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
# This script Scales the coefficients in a shortwave CKD look-up table
# in order to give the exact direct flux for a single reference
# atmosphere (usually the median atmosphere of the MMM dataset).

# Source the configuration and checking header scripts
. config.h
. check_configuration.h

SET=mmm
LBLINDIR=$MMM_SW_SPECTRA_DIR
LBLDIR=$WORK_SW_LBL_FLUX_DIR
LBLINPREFIX=$LBLINDIR/ckdmip_${SET}_sw_spectra_
CONFIG="config_sw_lbl_${SET}.nam"
BANDCODE=fluxes-raw
LBLPREFIX="ckdmip_${SET}_sw_${BANDCODE}"
SUFFIX=h5
ICOL=1
LBLFILE=$LBLDIR/${LBLPREFIX}_present_${ICOL}.$SUFFIX
SCENARIO=present

INDIR=${WORK_SW_RAW_CKD_DIR}
OUTDIR=${WORK_SW_RAW_CKD_DIR}

INCODE=raw-ckd-definition
OUTCODE=scaled-ckd-definition

# First perform a reference LBL calculation on the median/present gas
# concentrations, outputting direct fluxes only and omitting Rayleigh
# scattering.  Only do this if the file is not already present.
if [ ! -r "$LBLFILE" ]
then

  mkdir -p $LBLDIR

  cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
temperature_name = "temperature_hl",
do_write_spectral_fluxes = true,
do_write_optical_depth   = false,
do_write_direct_only     = true,
surf_albedo = 0.15,
use_mu0_dimension = false,
cos_solar_zenith_angle(1) = 0.5,
nblocksize = 1000,
iverbose = 3
/
EOF

  CO2_VMR=415e-6
  CH4_VMR=1921e-9
  N2O_VMR=332e-9

  H2O_FILE=${LBLINPREFIX}h2o_median.h5
  CO2_FILE=${LBLINPREFIX}co2_${SCENARIO}.h5
  O3_FILE=${LBLINPREFIX}o3_median.h5
  CH4_FILE=${LBLINPREFIX}ch4_${SCENARIO}.h5
  N2O_FILE=${LBLINPREFIX}n2o_${SCENARIO}.h5
  N2_FILE=${LBLINPREFIX}n2_constant.h5
  O2_FILE=${LBLINPREFIX}o2_constant.h5
  # RAYLEIGH_FILE=${LBLINPREFIX}rayleigh_${SCENARIO}_${COLS}.h5

  $CKDMIP_SW --config $CONFIG \
      --scenario "$SCENARIO" \
      --column-range $ICOL $ICOL \
      --ssi "$TRAINING_SW_SSI" \
      $H2O_FILE \
      $O3_FILE \
      $N2_FILE \
      $O2_FILE \
      --conc   $CO2_VMR   $CO2_FILE \
      --conc   $CH4_VMR   $CH4_FILE \
      --conc   $N2O_VMR   $N2O_FILE \
      $RAYLEIGH_FILE \
      --output $LBLFILE

fi


# Loop over each band structure and tolerance
for BANDSTRUCT in $BAND_STRUCTURE
do
    for TOL in $TOLERANCE
    do
	MODEL_CODE=${APPLICATION}_${BANDSTRUCT}-tol${TOL}${MODEL_CODE_SUFFIX}
	${BANNER} Scaling CKD model: $MODEL_CODE
	GPOINTFILE=${WORK_SW_GPOINTS_DIR}/${ECCKD_PREFIX}_sw_gpoints_${MODEL_CODE}.h5
	INPUT=${INDIR}/${ECCKD_PREFIX}_sw_${INCODE}_${MODEL_CODE}.nc
	OUTPUT=${OUTDIR}/${ECCKD_PREFIX}_sw_${OUTCODE}_${MODEL_CODE}.nc
	LOG=${OUTDIR}/${ECCKD_PREFIX}_sw_${OUTCODE}_${MODEL_CODE}.log

	$SCALE_LUT \
	    input=${INPUT} \
	    output=${OUTPUT} \
	    gpointfile=${GPOINTFILE} \
	    lblfile=${LBLFILE} \
	    model_id=sw_${APPLICATION}_${BANDSTRUCT}-tol${TOL} \
	    $EXTRA_ARGS \
	    |& tee $LOG
	test "${PIPESTATUS[0]}" -eq 0
    done
done
