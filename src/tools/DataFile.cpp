/// @file      DataFile.cpp
/// @brief     Implements the DataFile class
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#include <netcdf.h>

#include "DataFileEngineNetcdf.h"
#include "DataFileEngineCfg.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_GMVECSL
#include "DataFileEngineXml.h"
#endif

#include "DataFile.h"
#include "file_manager.h"

///
/// @details Firstly the absolute location of the file is found using
/// find_file(). If the filename starts with "." or "/" then it is
/// assumed to be relative to the current directory, or an absolute
/// path, respectively, in which case it is opened using
/// open_absolute(). Otherwise the list of data directories are
/// searched to find it. If the file is not found then
/// NO_PRODUCT_FOUND_ERROR is thrown.
///
void
DataFile::
open(const std::string& filename)
{
  std::string absfilename;
  absfilename = find_file(filename);
  if (absfilename.empty()) {
    ERROR << "Failed to find " << filename << " in " << search_directories();
    THROW(MISSING_MANDATORY_FILE);
  }
  else {
    open_absolute(absfilename);
  }
}

/// Close file if open
void
DataFile::
close()
{
  if (engine_) {
    delete engine_;
    engine_ = 0;
  }
}

///
/// @details Any existing file is first closed. If the filename has
/// extension ".nc", ".NC", ".cdf" or ".CDF" then the file is assumed
/// to be NetCDF format. If the code was compiled with NetCDF-4/HDF-5
/// support then the file extensions ".h5", ".hdf5" or ".HDF" then the
/// file is also assumed to be NetCDF format. A file extension of
/// ".cfg" is treated as CFG format. If the code was compiled against
/// the GMV Support Libraries then while an extension of ".xml",".XML" 
/// or ".EOF" is treated as XML format. An instance of the appropriate
/// DataFileEngine type is created, and subsequent "read" calls are
/// delegated to it. If the format is not recognised then a
/// PRODUCT_FORMAT_ERROR integer exception is thrown.
///
void
DataFile::
open_absolute(const std::string& filename)
{
  // Close existing file first
  close();

  // Determine type from filename
  size_t dotpos = filename.find_last_of('.');
  if (dotpos == std::string::npos) {
    ERROR << "Cannot determine file type for \"" << filename << "\"";
    THROW(PRODUCT_FORMAT_ERROR);
  }
  std::string extension = filename.substr(dotpos+1);
  if (extension == "nc" || extension == "NC"
      || extension == "cdf" || extension == "CDF"
#ifdef NC_NETCDF4
      || extension == "h5" || extension == "hdf5"
      || extension == "HDF"
#endif
      ) {
    engine_ = new DataFileEngineNetcdf(filename, throw_exceptions_);
  }
  else if (extension == "cfg") {
    engine_ = new DataFileEngineCfg(filename, throw_exceptions_);
  }
#ifdef HAVE_GMVECSL
  else if (extension == "xml" || extension == "XML" || extension == "EOF") {
    engine_ = new DataFileEngineXml(filename, throw_exceptions_);
  }
#endif
  else {
    ERROR << "File extension \"" << extension << "\" not recognised";
    THROW(PRODUCT_FORMAT_ERROR);
  }
}

///
/// @details This function takes as arguments those provided to the
/// main() function that describe the number and content of the
/// command-line arguments. The data file type is assumed to be CFG
/// and so it passes its arguments directly to the appropriate
/// constructor of the DataFileEngineCfg class.
///
bool
DataFile::
init(int argc, const char** argv) {
  close();
  engine_ = new DataFileEngineCfg(argc, argv, throw_exceptions_);
  return true;
}
