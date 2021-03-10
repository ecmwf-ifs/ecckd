#!/bin/bash

# To train CKD models for climate, we first optimize four gases: H2O,
# O3, CO2 and the composite of O2 and N2.  This is done by comparing
# to LBL calculations with the other gases set to zero.

# At ECMWF to have access to NCO tools
module load nco

. config.h

unset OMP_NUM_THREADS

SET=evaluation1

INDIR=$TRAINING_SW_SPECTRA_DIR
OUTDIR=$WORK_SW_LBL_FLUX_DIR

INPREFIX=$INDIR/ckdmip_${SET}_sw_spectra_

CONFIG="config_sw_lbl_evaluation.nam"

PROGRAM=$CKDMIP_SW

BANDCODE=fluxes
BANDCODE=fluxes-rgb

OUTPREFIX="ckdmip_${SET}_sw_${BANDCODE}"
SUFFIX=h5

STRIDE=1

cat > $CONFIG <<EOF
&shortwave_config
optical_depth_name = "optical_depth",
wavenumber_name = "wavenumber",
pressure_name = "pressure_hl",
temperature_name = "temperature_hl",
do_write_spectral_fluxes = false,
do_write_optical_depth   = false,
surf_albedo = 0.15,
use_mu0_dimension = true,
cos_solar_zenith_angle(1:5) = 0.1, 0.3, 0.5, 0.7, 0.9,
nspectralstride = $STRIDE,
nblocksize = 1000,
!band_wavenumber1(1:13) = 250, 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000,
!band_wavenumber2(1:13) = 2600, 3250, 4000, 4650, 5150, 6150, 8050, 12850, 16000, 22650, 29000, 38000, 50000,
!band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14286, 16667, 20000, 25000, 31746,
!band_wavenumber2(1:9) = 2500, 4000, 8000, 14286, 16667, 20000, 25000, 31746, 50000,
band_wavenumber1(1:9) = 250, 2500, 4000, 8000, 14300, 16650, 20000, 25000, 31750,
band_wavenumber2(1:9) = 2500, 4000, 8000, 14300, 16650, 20000, 25000, 31750, 50000,
iverbose = 3
/
EOF

mkdir -p $OUTDIR

SCENARIOS="rel-180 rel-280 rel-415 rel-560 rel-1120 rel-2240"

for SCENARIO in $SCENARIOS
do

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
	
	H2O_FILE=${INPREFIX}h2o_present_${COLS}.h5
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
	    $PROGRAM --config $CONFIG \
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
	    $PROGRAM --config $CONFIG \
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
	    $PROGRAM --config $CONFIG \
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
