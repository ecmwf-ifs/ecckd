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

unset OMP_NUM_THREADS

# The band code defines what bands the fluxes will be computed in,
# needed for training different band models
if [ "$#" -ge 1 ]
then
    BANDCODE=$1
else
    echo "Usage:"
    echo "  $0 <band>"
    echo 'where <band> can be fluxes (for the "narrow" band structure), fluxes-rgb, fluxes-fine or fluxes-vfine'
    exit 1
fi
#BANDCODE=fluxes-rgb
#BANDCODE=fluxes-fine
#BANDCODE=fluxes-vfine

SET=evaluation1

# This scenario only will also include spectral boundary fluxes
SPECTRAL_SCENARIO=rel-415

# Either use default water vapour continuum
H2OCONTINUUM=
# ...or another continuum model
#H2OCONTINUUM=mt-ckd-4.1.1

INDIR=${CKDMIP_DATA_DIR}/${SET}/sw_spectra

OUTDIR=$WORK_SW_LBL_FLUX_DIR

INPREFIX=$INDIR/ckdmip_${SET}_sw_spectra_

# H2O prefix to accommodate different continuum models
if [ ! "$H2OCONTINUUM" ]
then
    # Default continuum
    H2OPREFIX=$INDIR/ckdmip_${SET}_sw_spectra_
    H2OSUFFIX=
else
    H2OPREFIX=$WORK_SW_SPECTRA_DIR/ckdmip_${SET}_sw_spectra_
    H2OSUFFIX=-$H2OCONTINUUM
fi

CONFIG="config_sw_lbl_evaluation_pid$$.nam"
CONFIG_SPECTRAL="config_spectral_sw_lbl_evaluation_pid$$.nam"

PROGRAM=$CKDMIP_SW

OUTPREFIX="ckdmip_${SET}${H2OSUFFIX}_sw_${BANDCODE}"
SUFFIX=h5

STRIDE=1

# Specify bands corresponding to each band model
if [ "$BANDCODE" = fluxes ]
then
    BAND_WAVENUMBER1="band_wavenumber1(1:13) = 250, 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000,"
    BAND_WAVENUMBER2="band_wavenumber2(1:13) = 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000, 50000,"
elif [ "$BANDCODE" = fluxes-rgb ]
then
    BAND_WAVENUMBER1="band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14300, 16650, 20000, 25000, 31750,"
    BAND_WAVENUMBER2="band_wavenumber2(1:9) = 2500, 4000, 8000, 14300, 16650, 20000, 25000, 31750, 50000,"
elif [ "$BANDCODE" = fluxes-fine ]
then
    BAND_WAVENUMBER1="band_wavenumber1(1:26) = 250, 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 14300, 15400, 16650, 18200, 20000, 22200, 25000, 28550, 30250, 30750, 31250, 31750, 32250, 32750, 33250, 33750, 34250,"
    BAND_WAVENUMBER2="band_wavenumber2(1:26) = 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 14300, 15400, 16650, 18200, 20000, 22200, 25000, 28550, 30250, 30750, 31250, 31750, 32250, 32750, 33250, 33750, 34250, 50000,"
elif [ "$BANDCODE" = fluxes-vfine ]
then
    BAND_WAVENUMBER1="band_wavenumber1(1:44) = 250, 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 13800, 14300, 14800, 15400, 16000, 16650, 17400, 18200, 19050, 20000, 21050, 22200, 23550, 25000, 26300, 26650, 27050, 27400, 27800, 28150, 28550, 29000, 29400, 29850, 30300, 30750, 31250, 31750, 32250, 32800, 33350, 33900, 34500, 35100, 35700,"
    BAND_WAVENUMBER2="band_wavenumber2(1:44) = 2600, 3750, 5350, 7150, 8700, 10650, 12100, 13350, 13800, 14300, 14800, 15400, 16000, 16650, 17400, 18200, 19050, 20000, 21050, 22200, 23550, 25000, 26300, 26650, 27050, 27400, 27800, 28150, 28550, 29000, 29400, 29850, 30300, 30750, 31250, 31750, 32250, 32800, 33350, 33900, 34500, 35100, 35700, 50000,"
else
    echo "BANDCODE=$BANDCODE not understood"
    exit 1
fi

echo "Using BANDCODE=$BANDCODE"


# Standard namelist with band output
cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
temperature_name = "temperature_hl",
do_write_spectral_fluxes = false,
do_write_spectral_boundary_fluxes = false, ! Note !
do_write_optical_depth   = false,
surf_albedo = 0.15,
use_mu0_dimension = true,
cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
nspectralstride = $STRIDE,
nblocksize = 1000,
$BAND_WAVENUMBER1
$BAND_WAVENUMBER2
iverbose = 3
/
EOF

# Spectral namelist also outputting spectral boundary fluxes
cat > $CONFIG_SPECTRAL <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
temperature_name = "temperature_hl",
do_write_spectral_fluxes = false,
do_write_spectral_boundary_fluxes = true, ! Note !
do_write_optical_depth   = false,
surf_albedo = 0.15,
use_mu0_dimension = true,
cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
nspectralstride = $STRIDE,
nblocksize = 1000,
$BAND_WAVENUMBER1
$BAND_WAVENUMBER2
iverbose = 3
/
EOF

mkdir -p $OUTDIR

# ecCKD requires the following scenarios: for the first optimization
# step, well-mixed GHGs are set to CONSTANT except for CO2
SCENARIOS_REL="rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240"
# ...and in the subsequent steps the other gases are optimized
SCENARIOS_CKDMIP="present ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540"
# Combine the scenarios
SCENARIOS="$SCENARIOS_REL $SCENARIOS_CKDMIP"

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
# Ignore CFCs in the shortwave
#    CFC11_VMR=861e-12
#    CFC12_VMR=495e-12
    CFC11_VMR=0.0e-12
    CFC12_VMR=0.0e-12

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

	RAYLEIGH_FILE=${INPREFIX}rayleigh_present_${COLS}.h5

	OUTFILE=${OUTDIR}/RAW_${OUTPREFIXFULL}_${COLS}.${SUFFIX}
	OUTFILES="$OUTFILES $OUTFILE"
	
	if [ "$ONLY5GASES" = yes ]
	then
	    $PROGRAM --config $CONFIG_LOCAL \
		--scenario "$SCENARIO" \
		--ssi "$TRAINING_SW_SSI" \
		$H2O_FILE \
		$O3_FILE \
		$N2_ARG  $N2_FILE \
		$O2_ARG  $O2_FILE \
		--conc   $CO2_VMR $CO2_FILE \
		$RAYLEIGH_FILE \
		--output $OUTFILE
	elif [ "$RELATIVE" = yes ]
	then
	    $PROGRAM --config $CONFIG_LOCAL \
		--scenario "$SCENARIO" \
		--ssi "$TRAINING_SW_SSI" \
		$H2O_FILE \
		$O3_FILE \
		$N2_ARG  $N2_FILE \
		$O2_ARG  $O2_FILE \
		--const  $CH4_VMR   $CH4_FILE \
		--const  $N2O_VMR   $N2O_FILE \
		--conc   $CO2_VMR   $CO2_FILE \
		$RAYLEIGH_FILE \
		--output $OUTFILE
	else
	    $PROGRAM --config $CONFIG_LOCAL \
		--scenario "$SCENARIO" \
		--ssi "$TRAINING_SW_SSI" \
		$H2O_FILE \
		$O3_FILE \
		$N2_ARG  $N2_FILE \
		$O2_ARG  $O2_FILE \
		--conc   $CO2_VMR   $CO2_FILE \
		--conc   $CH4_VMR   $CH4_FILE \
		--conc   $N2O_VMR   $N2O_FILE \
		--conc   $CFC11_VMR $CFC11_FILE \
		--conc   $CFC12_VMR $CFC12_FILE \
		$RAYLEIGH_FILE \
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

rm -f $CONFIG $CONFIG_SPECTRAL
