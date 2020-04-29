# -*- shell-script -*-
# This script is included by the some of the others in this directory,
# and checks the inputs to those scripts. It also sets the APP and
# MIN_PRESSURE variables.

# Check necessary variables are present
if [ ! "${TOLERANCE}" ]
then
    ${BANNER_ERROR} '"TOLERANCE"' not specified
    exit 1
fi

if [ ! ${APPLICATION} ]
then
    ${BANNER_ERROR} '"APPLICATION"' not specified
    exit 1
fi

if [ ! "${BAND_STRUCTURE}" ]
then
    ${BANNER_ERROR} '"BAND_STRUCTURE"' not specified
    exit 1
fi

# Application-dependent actions
if [ "$APPLICATION" = climate ]
then 
    APP=climate
    MIN_PRESSURE=2
elif [ "$APPLICATION" = global-nwp ]
then
    APP=nwp
    MIN_PRESSURE=2
elif [ "$APPLICATION" = limited-area-nwp ]
then
    APP=nwp
    MIN_PRESSURE=400
else
    ${BANNER_ERROR} 'APPLICATION "'$APPLICATION'" not understood'
    exit 1
fi
