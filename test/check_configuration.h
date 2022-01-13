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
# This script is included by the some of the others in this directory,
# and checks the inputs to those scripts. It also sets the APP (if not
# already set) and MIN_PRESSURE variables.

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
    APP_LOCAL=climate
    MIN_PRESSURE=2
elif [ "$APPLICATION" = global-nwp ]
then
    APP_LOCAL=nwp
    MIN_PRESSURE=2
elif [ "$APPLICATION" = limited-area-nwp ]
then
    APP_LOCAL=nwp
    MIN_PRESSURE=400
else
    ${BANNER_ERROR} 'APPLICATION "'$APPLICATION'" not understood'
    exit 1
fi

# Only write APP if not already defined
if [ ! "$APP" ]
then
    APP=${APP_LOCAL}
fi
