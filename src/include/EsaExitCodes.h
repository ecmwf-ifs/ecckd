/// @file      EsaExitCodes.h
/// @brief     Program exit codes specified by ESA

#ifndef _EsaExitCodes_h
#define _EsaExitCodes_h

#include <map>
#include <string>
#include <iostream>

/* The exit codes are integers between 0 and 255 as defined in the
  "EarthCARE PDGS Generic IPF Interface Specifications" 
   EACA-GSEG-EOPG-TN-15-0001. The following exit codes have been pre-defined by
   ESA. */

#define EXIT_OK                         0
#define JOB_ORDER_READ_ERROR            128
#define PREMATURE_TERMINATION           129
#define OUT_OF_MEMORY                   130
#define UNEXPECTED_EXCEPTION            131
#define DISK_FULL                       132
#define XML_ERROR                       133
#define XML_WARNING                     1
#define MISSING_MANDATORY_FILE          134
#define NO_MATCHES_FOR_PATTERN          2
#define DUPLICATE_FILE                  3
#define INCONSISTENT_SPACECRAFT_ID      4
#define UNKNOWN_PROCESSING_TYPE         5
#define NOT_ENOUGH_INPUT_ERROR          135
#define NOT_ENOUGH_INPUT_WARNING        6
#define CRITICAL_FRAME_FAILURE          136
#define RECOVERABLE_FRAME_FAILURE       7
#define PRODUCT_FORMAT_ERROR            137
#define PRODUCT_FORMAT_WARNING          8
#define NO_PRODUCT_FOUND_ERROR          138
#define NO_PRODUCT_FOUND_WARNING        9
#define CANNOT_OPEN_MANDATORY_FILE      139
#define CANNOT_OPEN_OPTIONAL_FILE       10
#define NO_AUX_FILE_ERROR               140
#define NO_AUX_FILE_WARNING             11
#define BAD_MANDATORY_AUX_FORMAT        141
#define BAD_OPTIONAL_AUX_FORMAT         12
#define NO_PRODUCT_MODEL                142
#define READ_ERROR_PRODUCT_MODEL        143
#define FINAL_PRODUCT_CREATION_ERROR    144
#define TEMP_FILE_CREATION_ERROR        145
#define WRITE_ERROR                     146
#define PARAMETER_ERROR                 147
#define PARAMETER_WARNING               13
#define PROCESSING_ERROR                148
#define PROCESSING_WARNING              14

#define JOB_ORDER_READ_ERROR_STR         "Job order error"
#define PREMATURE_TERMINATION_STR        "Premature termination"
#define OUT_OF_MEMORY_STR                "Out of memory"
#define UNEXPECTED_EXCEPTION_STR         "Unexpected exception"
#define DISK_FULL_STR                    "Disk full"
#define XML_ERROR_STR                    "XML error"
#define XML_WARNING_STR                  "XML warning"
#define MISSING_MANDATORY_FILE_STR       "Missing mandatory file"
#define NO_MATCHES_FOR_PATTERN_STR       "No matches for pattern"
#define DUPLICATE_FILE_STR               "Duplicate file"
#define INCONSISTENT_SPACECRAFT_ID_STR   "Inconsistent spacecraft id"
#define UNKNOWN_PROCESSING_TYPE_STR      "Unknown processing type"
#define NOT_ENOUGH_INPUT_ERROR_STR       "Mandatory input data missing"
#define NOT_ENOUGH_INPUT_WARNING_STR     "Optional input data missing"
#define CRITICAL_FRAME_FAILURE_STR       "Critical frame failure"
#define RECOVERABLE_FRAME_FAILURE_STR    "Recoverable frame failure"
#define PRODUCT_FORMAT_ERROR_STR         "Product format error"
#define PRODUCT_FORMAT_WARNING_STR       "Product format warning"
#define NO_PRODUCT_FOUND_ERROR_STR       "No mandatory product found"
#define NO_PRODUCT_FOUND_WARNING_STR     "No optional product found"
#define CANNOT_OPEN_MANDATORY_FILE_STR   "Cannot open mandatory file"
#define CANNOT_OPEN_OPTIONAL_FILE_STR    "Cannot open optional file"
#define NO_AUX_FILE_ERROR_STR            "No mandatory aux file"
#define NO_AUX_FILE_WARNING_STR          "No optional aux file"
#define BAD_MANDATORY_AUX_FORMAT_STR     "Bad mandatory aux format"
#define BAD_OPTIONAL_AUX_FORMAT_STR      "Bad optional aux format"
#define NO_PRODUCT_MODEL_STR             "No product model"
#define READ_ERROR_PRODUCT_MODEL_STR     "Read error product model"
#define FINAL_PRODUCT_CREATION_ERROR_STR "Final product creation error"
#define TEMP_FILE_CREATION_ERROR_STR     "Temporary file creation error"
#define WRITE_ERROR_STR                  "Write error"
#define PARAMETER_ERROR_STR              "Parameter error"
#define PARAMETER_WARNING_STR            "Parameter warning"
#define PROCESSING_ERROR_STR             "Processing error"
#define PROCESSING_WARNING_STR           "Processing warning"

#endif  /* _EsaExitCodes_h */


