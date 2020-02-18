/// @file      OutputFileSettings.h
/// @brief     Preprocessor macros for output data
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright 2015-2016 European Centre for Medium Range Weather Forecasts
/// @license   Apache License Version 2 or ESA Software Community License Type 1
///            (see the NOTICE.md file for details)

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
