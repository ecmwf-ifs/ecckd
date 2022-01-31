/// @file      DataFile.h
/// @brief     Declares the class DataFile
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef DataFile_H
#define DataFile_H

#include "DataFileEngine.h"
#include "EsaExitCodes.h"

/// Used locally within this header file to check that a file has been
/// opened before trying to read from it
#define DATA_FILE_READ_CHECK if (!engine_) {				\
    ERROR << "Attempt to read from DataFile object not attached to a file"; \
    THROW(NO_PRODUCT_FOUND_ERROR); }


class OutputDataFile;

/// Provides the capability to read from configuration and data files
/// with different formats.
///
/// This class provides the ability to read strings, scalars, and
/// arrays of different numerical types from NetCDF, CFG (an ASCII
/// format) and XML files. Users of the class do not need to know
/// about the underlying format of the data. When a file is opened,
/// its type is determined from the file extension, and subsequent
/// access commands are implemented by delegating to the
/// DataFileEngine member object, an abstract base class from which
/// concrete classes inherit that implement the access commands for
/// each file type. Note that XML files may only be accessed if
/// compiling with the GMV Support Library.
///
/// Many functions have two forms, one in which only the variable name
/// is provided, and another in which a "scope" is specified as well
/// as the "varname". In the case of CFG files, the scope represents a
/// named section of the file.  For NetCDF files this form is a way to
/// extract attributes: "scope" is interpreted as a NetCDF variable
/// name and "varname" as the attribute name. Use
/// DATA_FILE_GLOBAL_SCOPE as the scope argument to get the global
/// attributes for a NetCDF file.
class DataFile
{
public:
  /// Construct an object unattached to a file
 DataFile() : engine_(0), throw_exceptions_(false) { }

  /// Constructor that immediately opens a file
  DataFile(std::string filename, ///< Name of file to open
	   bool throw_exceptions = false) ///< Do we throw exceptions on read failure of just return "false"
    : engine_(0), throw_exceptions_(throw_exceptions)
  { open(filename); }

  /// Constructor that reads a CFG file from the command line
  /// arguments and supplements with any param=value pairs
  DataFile(int argc, const char** argv)
    : engine_(0), throw_exceptions_(false)
  { init(argc, argv); }

  /// Destructor closes the file if open
  ~DataFile() { close(); }

  /// @name Opening and closing files
  /// @{

  /// Open a data file for reading using the search path to find it,
  /// unless the filename begins with "." or "/"
  void open(const std::string& filename);

  /// As "open" but with an absolute file name
  void open_absolute(const std::string& filename);

  /// Read a CFG file from the command line arguments and supplement
  /// with any param=value pairs
  bool init(int argc, const char** argv);

  /// Close file if open
  void close();

  /// @}

  /// @name Inquiry functions
  /// @{

  /// Return true if a file is open, false otherwise
  bool is_open() const { return (engine_ != 0); }

  /// Return the name of the open file
  std::string name() const {
    DATA_FILE_READ_CHECK;
    return engine_->name();
  }

  /// Return the length of the dimensions of a variable "varname"
  IntVector size(const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->size(varname);
  }

  /// Return the dimensions of a variable "varname" in the specified
  /// "scope"
  IntVector size(const std::string& scope, const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->size(scope, varname);
  }

  /// Return "true" if a variable exists, "false" otherwise
  bool exist(const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->exist(varname);
  }

  /// Return "true" if the variable "varname" exists in "scope"
  bool exist(const std::string& scope, const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->exist(scope, varname);
  }

  /// Get a list of variable or argument names present in the file,
  /// where global scope is assumed if "scope" is omitted
  bool get_names(std::list<std::string>& names,
		 const std::string& scope = std::string()) const {
    DATA_FILE_READ_CHECK;
    return engine_->get_names(names, scope);
  }

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
  bool read(Real& x, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(x, varname, j, i);
  }

  /// Read a real scalar from the specified "scope" of the data file
  bool read(Real& x, const std::string& scope,
	    const std::string& varname, int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(x, scope, varname, j, i);
  }

  /// Read an integer scalar from data file
  bool read(int& x, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(x, varname, j, i);
  }

  /// Read an integer scalar from the specified "scope" of the data file
  bool read(int& x, const std::string& scope,
	    const std::string& varname, int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(x, scope, varname, j, i);
  }

  /// Read a boolean from the data file
  bool read(bool& x, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(x, varname, j, i);
  }

  /// Read a boolean scalar from the specified "scope" of the data file
  bool read(bool& x, const std::string& scope,
	    const std::string& varname, int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(x, scope, varname, j, i);
  }

  /// Read a string from the data file
  bool read(std::string& s, const std::string& varname,
	    int j = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(s, varname, j);
  }

  /// Read a string from the specified "scope" of the data file
  bool read(std::string& s, const std::string& scope,
	    const std::string& varname, int j = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(s, scope, varname, j);
  }

  /// Read a vector of real numbers from the data file
  bool read(Vector& v, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(v, varname, j, i);
  }

  /// Read a slice real numbers from the data file
  bool read_slice(Vector& v, const std::string& varname,
	    int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read_slice(v, varname, i);
  }

  /// Read a real vector from the specified "scope" of the data file
  bool read(Vector& v, const std::string& scope,
	    const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(v, scope, varname);
  }

  /// Read a vector of integers from the data file
  bool read(IntVector& v, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(v, varname, j, i);
  }

  /// Read a vector of integers from the specified "scope" of the data
  /// file
  bool read(IntVector& v, const std::string& scope,
	    const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(v, scope, varname);
  }

  /// Read a real matrix from the data file
  bool read(Matrix& M, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(M, varname, j, i);
  }

  /// Read a real matrix from the specified "scope of the data file
  bool read(Matrix& M, const std::string& scope,
	    const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(M, scope, varname);
  }

  /// Read a real 3D array from the data file
  bool read(Array3D& A, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(A, varname, j, i);
  }
  /*
  bool read(Array3& A, const std::string& scope,
	    const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(A, scope, varname);
  }
  */

  /// Read a real 4D array from the data file
  bool read(Array4D& A, const std::string& varname,
	    int j = -1, int i = -1) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(A, varname, j, i);
  }
  /// Read the entire contents of the file to a string
  bool read(std::string& s) const {
    DATA_FILE_READ_CHECK;
    return engine_->read(s);
  }

  /// Read the missing_value associated with a varname
  bool read_missing_value(Real& missing_val, const std::string& varname) const {
    DATA_FILE_READ_CHECK;
    return engine_->read_missing_value(missing_val, varname);
  }

  ///@}

  /// @name Miscellaneous
  /// @{

  /// If the argument is true then any failed read will throw an
  /// exception
  void throw_exceptions(bool in) { throw_exceptions_ = in; }

  /// Copy all the attributes of varname to a variable in the specified
  /// output file; if the last argument is omitted it will be assumed
  /// that the target variable has the same name as the source variable
  bool copy_attributes(const std::string& varname,
		       OutputDataFile* output,
		       const std::string& output_varname = std::string()) 
    const {
    DATA_FILE_READ_CHECK;
    return engine_->copy_attributes(varname, output, output_varname);
  }

  /// Write the configuration variables that have been read and logged
  /// to the specified file, returning true on success, false
  /// otherwise (e.g. if the input was not capable of logging the
  /// variables read)
  bool write_stored_variables(const std::string& filename) const {
    DATA_FILE_READ_CHECK;
    return engine_->write_stored_variables(filename);
  }

  /// @}

protected:
  /// Pointer to the object that provides access to the particular
  /// format of data file
  DataFileEngine* engine_;

  /// Do we throw exceptions if an error occurs reading a variable?
  bool throw_exceptions_;

};

#undef DATA_FILE_READ_CHECK

#endif
