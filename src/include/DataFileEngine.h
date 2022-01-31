/// @file      DataFileEngine.h
/// @brief     Declares the class DataFileEngine
/// @author    Robin J. Hogan
/// @copyright (C) Copyright 2016- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef DataFileEngine_H
#define DataFileEngine_H

#include <string>
#include <list>

#include "arrays.h"
/*
enum DataFileType {
  UNDETERMINED,
  CONFIG, 
  NETCDF
};
*/

/// Use this as the "scope" argument of the read commands to access
/// the global attributes of a NetCDF file
#define DATA_FILE_GLOBAL_SCOPE "_global_"

class OutputDataFile;

/// Base class providing the template for reading different types of
/// data file
///
/// This class is used by the DataFile class.  Many of its functions
/// have two forms, one in which only the variable name is provided,
/// and another in which a "scope" is specified as well as the
/// "varname". In the case of CFG files, the scope represents a named
/// section of the file.  For NetCDF files this form is a way to
/// extract attributes: "scope" is interpreted as a NetCDF variable
/// name and "varname" as the attribute name.
class DataFileEngine
{
public:

  /// Close the file
  virtual ~DataFileEngine() {};

  /// @name Inquiry functions
  /// @{

  /// Return the name of the open file
  std::string name() const { return filename_; }

  /// Return the length of the dimensions of a variable "varname"
  virtual IntVector size(const std::string& varname) const = 0;

  /// Return the dimensions of a variable "varname" in the specified
  /// "scope"
  virtual IntVector size(const std::string& scope, const std::string& varname) const = 0;

  /// Return "true" if a variable exists, "false" otherwise
  virtual bool exist(const std::string& varname) const = 0;

  /// Return "true" if the variable "varname" exists in "scope"
  virtual bool exist(const std::string& scope, const std::string& varname) const = 0;

  /// Get a list of variable or argument names present in the file,
  /// where global scope is assumed if "scope" is omitted
  virtual bool get_names(std::list<std::string>& names,
			 const std::string& scope = std::string()) const = 0;

  /// @}

  /// @name Reading data
  ///
  /// @details In the following functions, if "j" and "i", or just
  /// "i", are provided then a slice or element is being requested
  /// from the assumed larger dimensional array in the data file. For
  /// example, if a scalar is requested and just "j" is provided then
  /// the variable is assumed to be a vector and element "j" is being
  /// requested (0-based indexing). If a scalar is requested and both
  /// "j" and "i" are provided then the variable in the file is
  /// assumed to be a matrix and the element at row j, column i is
  /// requested. If a vector is requested and just "j" is provided
  /// then the array in the file is assumed to be a matrix and slice
  /// "j" is requested, with j indexing the slowest-varying dimension
  /// in memory.
  ///
  /// If "throw_exceptions" is false then these functions return
  /// "true" if successful, "false" otherwise (in which case x will
  /// not be touched). If "throw_exceptions" is true then the
  /// functions will throw an error if the requested variable was not
  /// found or if any other problem occurred.
  ///
  ///  @{

  /// Read a real scalar from the file, returning true if successful
  /// and false if not
  virtual bool read(Real& x, const std::string& varname,
		    int j = -1, int i = -1) const = 0;

  /// Read a real scalar from the specified "scope" of the data file
  virtual bool read(Real& x, const std::string& scope,
		    const std::string& varname, int j = -1, int i = -1) const = 0;
		  
  /// Read an integer scalar from data file
  virtual bool read(int& x, const std::string& varname,
		    int j = -1, int i = -1) const = 0;

  /// Read an integer scalar from the specified "scope" of the data file
  virtual bool read(int& x, const std::string& scope,
		    const std::string& varname, int j = -1, int i = -1) const = 0;

  /// Read a boolean from the data file
  virtual bool read(bool& x, const std::string& varname,
		    int j = -1, int i = -1) const = 0;

  /// Read a boolean scalar from the specified "scope" of the data file
  virtual bool read(bool& x, const std::string& scope,
		    const std::string& varname, int j = -1, int i = -1) const = 0;

  /// Read a string from the data file
  virtual bool read(std::string& s, const std::string& varname,
		    int j = -1) const = 0;

  /// Read a string from the specified "scope" of the data file
  virtual bool read(std::string& s, const std::string& scope,
		    const std::string& varname, int j = -1) const = 0;

  /// Read a slice of real numbers from the data file
  virtual bool read_slice(Vector& v, const std::string& varname,
		    int i = -1) const = 0;

  /// Read a vector of real numbers from the data file
  virtual bool read(Vector& v, const std::string& varname,
		    int j = -1, int i = -1) const = 0;

  /// Read a real vector from the specified "scope" of the data file
  virtual bool read(Vector& v, const std::string& scope,
		    const std::string& varname) const = 0;

  /// Read a vector of integers from the data file
  virtual bool read(IntVector& v, const std::string& varname,
		    int j = -1, int i = -1) const = 0;

  /// Read a vector of integers from the specified "scope" of the data
  /// file
  virtual bool read(IntVector& v, const std::string& scope,
		    const std::string& varname) const = 0;

  /// Read a real matrix from the data file
  virtual bool read(Matrix& M, const std::string& varname,
		 int j = -1, int i = -1) const = 0;

  /// Read a real matrix from the specified "scope of the data file
  virtual bool read(Matrix& M, const std::string& scope,
  		    const std::string& varname) const = 0;

  /// Read a real 3D array from the data file
  virtual bool read(Array3D& A, const std::string& varname,
		    int j = -1, int i = -1) const = 0;
  //  virtual bool read(Array3& A, const std::string& scope,
  //		    const std::string& varname) const = 0;

  virtual bool read(Array4D& A, const std::string& varname,
		    int j = -1, int i = -1) const = 0;

  /// Read the entire contents of the file to a string
  virtual bool read(std::string& s) const = 0;

  /// Read the missing_value associated with a varname
  virtual bool read_missing_value(Real& missing_val, const std::string& varname) const = 0;
	
  /// @}

  /// @name Miscellaneous
  /// @{

  /// Copy all the attributes of varname to a variable in the specified
  /// output file; if the last argument is omitted it will be assumed
  /// that the target variable has the same name as the source variable
  virtual bool copy_attributes(const std::string& varname,
			       OutputDataFile* output,
			       const std::string& output_varname = std::string())  
    const = 0;

  /// Write the configuration variables that have been read and logged
  /// to the specified file, returning true on success, false
  /// otherwise
  virtual bool write_stored_variables(const std::string& filename) const { return false; }

  /// @}

 protected:
  /// Name of the data file name
  std::string filename_;
};


#endif
