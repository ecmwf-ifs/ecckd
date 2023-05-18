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
# To train CKD models for climate, we first optimize four gases: H2O,
# O3, CO2 and the composite of O2 and N2.  This is done by comparing
# to LBL calculations with the other gases set to zero.

# At ECMWF to have access to NCO tools
module load nco

. config.h

SET=evaluation1

# This scenario only will also include spectral boundary fluxes
SPECTRAL_SCENARIO=rel-415

# Either use default water vapour continuum
H2OCONTINUUM=
# ...or another continuum model
#H2OCONTINUUM=mt-ckd-4.1.1

INDIR=${CKDMIP_DATA_DIR}/${SET}/lw_spectra

OUTDIR=$WORK_LW_LBL_FLUX_DIR

INPREFIX=$INDIR/ckdmip_${SET}_lw_spectra_


# H2O prefix to accommodate different continuum models
if [ ! "$H2OCONTINUUM" ]
then
    # Default continuum
    H2OPREFIX=$INDIR/ckdmip_${SET}_lw_spectra_
    H2OSUFFIX=
else
    H2OPREFIX=$WORK_LW_SPECTRA_DIR/ckdmip_${SET}_lw_spectra_
    H2OSUFFIX=-$H2OCONTINUUM
fi

CONFIG="config_lw_lbl_evaluation_pid$$.nam"
CONFIG_SPECTRAL="config_spectral_lw_lbl_evaluation_pid$$.nam"

PROGRAM=$CKDMIP_LW

OUTPREFIX="ckdmip_${SET}${H2OSUFFIX}_lw_fluxes"
SUFFIX=h5

STRIDE=1

# Number of angles per hemisphere (0=classic two-stream)
NANGLE=0

# Standard namelist with band output
cat > $CONFIG <<EOF
&longwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
nspectralstride = $STRIDE,
nangle          = $NANGLE,
do_write_planck = .false.,
do_write_spectral_fluxes = .false.,
do_write_spectral_boundary_fluxes = false, ! Note !
do_write_optical_depth = .false.,
band_wavenumber1(1:13) = 0, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080,
band_wavenumber2(1:13) = 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 3260,
iverbose = 3
/
EOF

# Spectral namelist also outputting spectral boundary fluxes
cat > $CONFIG_SPECTRAL <<EOF
&longwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
pressure_scaling = 1.0,
temperature_name = "temperature_hl",
nspectralstride = $STRIDE,
nangle          = $NANGLE,
do_write_planck = .false.,
do_write_spectral_fluxes = .false.,
do_write_spectral_boundary_fluxes = true, ! Note !
do_write_optical_depth = .false.,
band_wavenumber1(1:13) = 0, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080,
band_wavenumber2(1:13) = 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 3260,
iverbose = 3
/
EOF

# File prefix contains _fluxes_ for classic two-stream,
# _fluxes-Nangle_ for N angles per hemisphere
if [ "$NANGLE" -gt 0 ]
then
    OUTPREFIX=${OUTPREFIX}-${NANGLE}angle
fi

mkdir -p $OUTDIR

# ecCKD requires the following scenarios: for the first optimization
# step, well-mixed GHGs are set to CONSTANT except for CO2
SCENARIOS_REL="rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240"
# ...and in the subsequent steps the other gases are optimized
SCENARIOS_CKDMIP="present ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540"
SCENARIOS_CKDMIP_CFC="cfc11-0 cfc11-2000 cfc12-0 cfc12-550"
# Combine the scenarios
SCENARIOS="$SCENARIOS_REL $SCENARIOS_CKDMIP $SCENARIOS_CKDMIP_CFC"

for SCENARIO in $SCENARIOS
do
    CONFIG_LOCAL=$CONFIG
    # Is this the spectral scenario?  If so, use a different namelist
    # file
    if [ "$SCENARIO" = "$SPECTRAL_SCENARIO" ]
    then
	CONFIG_LOCAL=$CONFIG_SPECTRAL
    fi
    
    OUTPREFIXFULL=${OUTPREFIX}_${SCENARIO}

    # Default concentrations are present-day ones
    CO2_VMR=415e-6
    CH4_VMR=1921e-9
    N2O_VMR=332e-9
    CFC11_VMR=861e-12
    CFC12_VMR=495e-12

    O2_ARG=
    N2_ARG=

    ONLY5GASES=no
    RELATIVE=no

    NF=$(echo $SCENARIO | awk -F- '{print NF}')
    # Set concentrations of trace gases according to SCENARIO
    if [ "$NF" = 4 ]
    then
	# We have a combination of concentration perturbations
	GAS1=$(echo $SCENARIO | awk -F- '{print $1}')
	CONC1=$(echo $SCENARIO | awk -F- '{print $2}')
	GAS2=$(echo $SCENARIO | awk -F- '{print $3}')
	CONC2=$(echo $SCENARIO | awk -F- '{print $4}')
	if [ "$GAS1" = "co2" ]
	then
	    CO2_VMR=${CONC1}e-6
	elif [ "$GAS1" = "ch4" ]
	then
	    CH4_VMR=${CONC1}e-9
	elif [ "$GAS1" = "n2o" ]
	then
	    N2O_VMR=${CONC1}e-9
	else
	    echo "First gas \"$GAS1\" in combination not recognised"
	    exit
	fi
	if [ "$GAS2" = "co2" ]
	then
	    CO2_VMR=${CONC2}e-6
	elif [ "$GAS2" = "ch4" ]
	then
	    CH4_VMR=${CONC2}e-9
	elif [ "$GAS2" = "n2o" ]
	then
	    N2O_VMR=${CONC2}e-9
	else
	    echo "Second gas \"$GAS2\" in combination not recognised"
	    exit
	fi
    elif [ "$SCENARIO" = present ]
    then
	CO2_VMR=415e-6
	CH4_VMR=1921e-9
	N2O_VMR=332e-9
	CFC11_VMR=861e-12
	CFC12_VMR=495e-12
    elif [ "$SCENARIO" = preindustrial ]
    then
	CO2_VMR=280e-6
	CH4_VMR=700e-9
	N2O_VMR=270e-9
	CFC11_VMR=32e-12
	CFC12_VMR=0e-12
    elif [ "$SCENARIO" = future ]
    then
	CO2_VMR=1120e-6
	CH4_VMR=3500e-9
	N2O_VMR=405e-9
	CFC11_VMR=2000e-12
	CFC12_VMR=200e-12
    elif [ "$SCENARIO" = glacialmax ]
    then
	CO2_VMR=180e-6
	CH4_VMR=350e-9
	N2O_VMR=190e-9
	CFC11_VMR=32e-12
	CFC12_VMR=0e-12
    elif [ "${SCENARIO:0:4}" = co2- ]
    then
	CO2_VMR=${SCENARIO:4:100}e-6
    elif [ "${SCENARIO:0:5}" = 5gas- ]
    then
	# Only O2, N2, H2O, O3 and CO2
	CO2_VMR=${SCENARIO:5:100}e-6
	CH4_VMR="0.0"
	N2O_VMR="0.0"
	CFC11_VMR="0.0"
	CFC12_VMR="0.0"
	ONLY5GASES=yes
    elif [ "${SCENARIO:0:4}" = rel- ]
    then
	# CO2 allowed to vary; N2O and CH4 use constant mole fractions
	# as they are to be parameterized using the "relative-linear"
	# scheme
	CO2_VMR=${SCENARIO:4:100}e-6
	RELATIVE=yes
    elif [ "${SCENARIO:0:4}" = ch4- ]
    then
	CH4_VMR=${SCENARIO:4:100}e-9
    elif [ "${SCENARIO:0:4}" = n2o- ]
    then
	N2O_VMR=${SCENARIO:4:100}e-9
    elif [ "${SCENARIO:0:6}" = cfc11- ]
    then
	CFC11_VMR=${SCENARIO:6:100}e-12
    elif [ "${SCENARIO:0:6}" = cfc12- ]
    then
	CFC12_VMR=${SCENARIO:6:100}e-12
    elif [ "$SCENARIO" = o2-0 ]
    then
	O2_ARG="--scale 0"
    elif [ "$SCENARIO" = n2-0 ]
    then
	N2_ARG="--scale 0"
    else    
	echo "SCENARIO=$SCENARIO not understood"
	exit
    fi
    
    OUTFILES=""

    # Loop over 50 columns in groups of 10
    for STARTCOL in 1 11 21 31 41
    do
	ENDCOL=$(expr $STARTCOL + 9)
	COLS=${STARTCOL}-${ENDCOL}
	
	H2O_FILE=${H2OPREFIX}h2o${H2OSUFFIX}_present_${COLS}.h5
	CO2_FILE=${INPREFIX}co2_present_${COLS}.h5
	O3_FILE=${INPREFIX}o3_present_${COLS}.h5
	CH4_FILE=${INPREFIX}ch4_present_${COLS}.h5
	N2O_FILE=${INPREFIX}n2o_present_${COLS}.h5
	CFC11_FILE=${INPREFIX}cfc11_present-equivalent_${COLS}.h5
	CFC12_FILE=${INPREFIX}cfc12_present_${COLS}.h5
	N2_FILE=${INPREFIX}n2_constant_${COLS}.h5
	O2_FILE=${INPREFIX}o2_constant_${COLS}.h5
	
	OUTFILE=${OUTDIR}/RAW_${OUTPREFIXFULL}_${COLS}.${SUFFIX}
	OUTFILES="$OUTFILES $OUTFILE"
	
	if [ "$ONLY5GASES" = yes ]
	then
	    $PROGRAM --config $CONFIG_LOCAL \
		--scenario "$SCENARIO" \
		$H2O_FILE \
		$O3_FILE \
		$N2_ARG  $N2_FILE \
		$O2_ARG  $O2_FILE \
		--conc   $CO2_VMR $CO2_FILE \
		--output $OUTFILE
	elif [ "$RELATIVE" = yes ]
	then
	    $PROGRAM --config $CONFIG_LOCAL \
		--scenario "$SCENARIO" \
		$H2O_FILE \
		$O3_FILE \
		$N2_ARG  $N2_FILE \
		$O2_ARG  $O2_FILE \
		--const  $CH4_VMR   $CH4_FILE \
		--const  $N2O_VMR   $N2O_FILE \
		--conc   $CO2_VMR   $CO2_FILE \
		--output $OUTFILE
	else
	    $PROGRAM --config $CONFIG_LOCAL \
		--scenario "$SCENARIO" \
		$H2O_FILE \
		$O3_FILE \
		$N2_ARG  $N2_FILE \
		$O2_ARG  $O2_FILE \
		--conc   $CO2_VMR   $CO2_FILE \
		--conc   $CH4_VMR   $CH4_FILE \
		--conc   $N2O_VMR   $N2O_FILE \
		--conc   $CFC11_VMR $CFC11_FILE \
		--conc   $CFC12_VMR $CFC12_FILE \
		--output $OUTFILE
	fi
    done
    
    OUTFILE=${OUTDIR}/${OUTPREFIXFULL}.${SUFFIX}
    
    echo "*** WRITING $OUTFILE ***"

    ncrcat -O $OUTFILES $OUTFILE
    
    if [ "$?" = 0 ]
    then
	rm -f $OUTFILES
    fi

done
