/// @file      OutputDataFile.h
/// @brief     Declares the class OutputDataFile
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef OutputDataFile_H
#define OutputDataFile_H

#include <string>
#include <cstdio>

#include <adept_arrays.h>

#include "file_manager.h"
#include "readconfig.h"
#include "DataFile.h"
#include "OutputFileSettings.h"

enum OutputDataFileMode {
  OUTPUT_MODE_CLOBBER,
  OUTPUT_MODE_NOCLOBBER,
  OUTPUT_MODE_APPEND
};

enum OutputDataFileVariableType {
  DOUBLE,
  FLOAT,
  INT,
  SHORT,
  STRING,
  BYTE  
};

enum DataFileType {
  UNDETERMINED,
  CONFIG, 
  NETCDF
};

/// Provides the capability to write NetCDF and CFG data files
class OutputDataFile
{
public:
  //Constructor
  OutputDataFile() : type_(UNDETERMINED), is_netcdf4hdf5_(false), file_(0), ncid_(-1),
    throw_exceptions_(false),chunk_size_(-1),write_full_err_stats_(1) { }

  OutputDataFile(std::string filename, OutputDataFileMode mode = OUTPUT_MODE_CLOBBER,
		 bool throw_exceptions = false)
    : type_(UNDETERMINED),is_netcdf4hdf5_(false), file_(0), ncid_(-1),
    throw_exceptions_(throw_exceptions),chunk_size_(-1),write_full_err_stats_(1) 
  { open(filename, mode); }

  ~OutputDataFile() { close(); }

  // Open a data file for writing - note that it will use the first
  // search directory from file_manager.h
  void open(const std::string& filename, 
	    OutputDataFileMode = OUTPUT_MODE_CLOBBER,
	    DataFileType = UNDETERMINED);
  // As "open", but with an absolute file name
  void open_absolute(const std::string& filename, 
		     OutputDataFileMode = OUTPUT_MODE_CLOBBER,
		     DataFileType = UNDETERMINED);
  void close();

  void set_options(DataFile& config);

  // NetCDF files require the dimensions to be defined; if the length
  // is omitted then the dimension is assumed to be "unlimited"
  void define_dimension(const std::string& dimname, int length = -1);
  void define_variable(const std::string& varname, 
		       OutputDataFileVariableType type,
		       const std::string& dim1name = std::string(),
		       const std::string& dim2name = std::string(),
		       const std::string& dim3name = std::string(),
		       const std::string& dim4name = std::string());
  void deflate_variable(const std::string& varname);
  void write_units(const std::string& units, 
		   const std::string& varname) {
    write(units, varname, "units");
  }
  void write_long_name(const std::string& long_name,
		       const std::string& varname) {
    write(long_name, varname, "long_name");
  }
  void write_comment(const std::string& comment, 
		     const std::string& varname) {
    write(comment, varname, "comment");
  }
  bool write_missing_value(Real missing_value,
			   const std::string& varname);
  void write_plot_range(Real min, Real max, bool is_logarithmic,
			const std::string& varname) {
    Vector range(2);
    range(0) = min;
    range(1) = max;
    write(range, varname, "plot_range");
    if (is_logarithmic) {
      write("logarithmic", varname, "plot_scale");
    }
    else {
      write("linear", varname, "plot_scale");
    }
  }

  // Append the "title" attribute
  // (NetCDF only)
  bool append_title();
  
  // Append the "conventions" attribute
  // (NetCDF only)
  bool append_convention();

  // Append the current command line to the "history" attribute
  // (NetCDF only)
  bool append_history(int argc, const char* argv[]);
  
  void end_define_mode();

  bool write_comment(const std::string& comment);

  // Write a real number to data file
  bool write(Real x, const std::string& varname,
	    int j = -1, int i = -1);
  // Write a real number to specified "scope" of data file, which
  // could be a section in a cfg file or an attribute of a variable in
  // the case of a netcdf file
  bool write(Real x, const std::string& scope,
	    const std::string& varname, int j = -1, int i = -1);

  // Write an integer to data file
  bool write(int x, const std::string& varname,
	    int j = -1, int i = -1);
  bool write(int x, const std::string& scope,
	    const std::string& varname, int j = -1, int i = -1);

  // Write a string to a data file
  bool write(const std::string& s, const std::string& varname,
	    int j = -1);
  bool write(const std::string& s, const std::string& scope,
	    const std::string& varname, int j = -1);

  // Write a vector of real numbers to a data file
  bool write(const Vector& v, const std::string& varname,
	    int j = -1, int i = -1);
  bool write(const Vector& v, const std::string& scope,
	    const std::string& varname);

  // Write a vector of integers to a data file
  bool write(const IntVector& v, const std::string& varname,
	    int j = -1, int i = -1);
  bool write(const IntVector& v, const std::string& scope,
	    const std::string& varname);

  // Write a matrix of numbers to a data file
  bool write(const Matrix& M, const std::string& varname,
	    int j = -1, int i = -1);
  bool write(const Matrix& M, const std::string& scope,
	    const std::string& varname);

  // Write a 3D array of numbers to a data file
  bool write(const Array3D& M, const std::string& varname,
	    int j = -1, int i = -1);

  bool exist(const std::string& varname);
  bool exist(const std::string& scope, const std::string& varname);

  bool use_chunking() { return chunk_size_ > 0; }
  void set_chunking(const IntVector& chunk_size, const std::string& varname);

  bool write_full_error_stats() { return write_full_err_stats_ > 0; }

private:
  DataFileType type_;
  bool is_netcdf4hdf5_;
  FILE* file_;
  std::string filename_;
  std::string time_dimension_name_;
  std::string range_dimension_name_;
  int root_ncid_;
  int header_ncid_;
  int ncid_;
  bool throw_exceptions_;
  int chunk_size_;
  int write_full_err_stats_;
};



#endif
