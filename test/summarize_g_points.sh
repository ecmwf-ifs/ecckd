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
# This script summarizes the number of g points in the NetCDF files
# provided on the command line

for FILE in $@
do
    INFO=$(ncdump -h $FILE | sed -n '3,7p')
    echo $INFO $FILE
done
