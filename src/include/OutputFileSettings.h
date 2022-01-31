/// @file      OutputFileSettings.h
/// @brief     Preprocessor macros for output data
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef OutputFileSettings_H
#define OutputFileSettings_H 1

#include <netcdf.h>

//#define OUTPUT_TIME_DIMENSION "time"
extern std::string OUTPUT_TIME_DIMENSION;
#define OUTPUT_RANGE_DIMENSION_NAME "height"
extern std::string OUTPUT_RANGE_DIMENSION;
#define OUTPUT_ITERATION_DIMENSION "iteration"
#define OUTPUT_THERMODYNAMIC_RANGE_DIMENSION "thermodynamic_height"

#define OUTPUT_MISSING_VALUE NC_FILL_FLOAT
//-999.0

#define OUTPUT_REAL_TYPE FLOAT

#define OUTPUT_COST_FUNCTION_NAME "cost_function"

#define OUTPUT_OBSERVATIONS 1

#define OUTPUT_TRUE_SUFFIX "_true"
#define OUTPUT_ITERATION_SUFFIX "_iteration"

#endif
