/// @file      DataFileEngineXml.h
/// @brief     Declares the class DataFileEngineXml
/// @author    Robin J. Hogan
/// @copyright (C) Copyright 2016- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef DataFileEngineXml_H
#define DataFileEngineXml_H

#include <map>
#include <string>

#include "DataFileEngine.h"

/// Provides a containter to store a variable when it is read
class DataLogXml
{
public:
  void save(const char* s) { str_ = s; type_ = LOG_STRING; }
  void save(const std::string& s) { str_ = s; type_ = LOG_STRING; }
  void save(const Vector& v) {
    vector_.resize(v.size());
    for (int i = 0; i < v.size(); ++i) { vector_[i] = v[i]; }
    type_ = LOG_VECTOR;
  }
  void save(const std::vector<double>& v) { vector_ = v; type_ = LOG_VECTOR; }
  void save(const IntVector& v) {
    int_vector_.resize(v.size());
    for (int i = 0; i < v.size(); ++i) { int_vector_[i] = v[i]; }
    type_ = LOG_INT_VECTOR;
  }
  void save(const std::vector<long>& v) { int_vector_ = v; type_ = LOG_INT_VECTOR; }
  void save(double x) { real_ = x; type_ = LOG_REAL; }
  void save(long i) { int_ = i; type_ = LOG_INT; }

  /// Write this variable to an XML output file, using the appropriate type
  int write(const std::string& filename, const std::string& varname) const;

private:
  enum { LOG_STRING, LOG_VECTOR, LOG_INT_VECTOR,
	 LOG_REAL, LOG_INT } type_;
  std::string str_;
  std::vector<double> vector_;
  std::vector<long> int_vector_;
  double real_;
  long int_;
};


/// Provides the capability for users of the DataFile class to read
/// XML files
///
/// Note that this capability is only provided if the code is compiled
/// against the GMV Support Library
class DataFileEngineXml : public DataFileEngine
{
public:

  /// Open a CFG file with an absolute path and specify whether read
  /// errors should throw an exception or just return "false".
  DataFileEngineXml(const std::string& filename,
		    bool throw_exceptions = false);

  /// Close the file
  virtual ~DataFileEngineXml();

  /// @name Inquiry functions
  /// @{

  /// Return the length of the dimensions of a variable "varname"
  virtual IntVector size(const std::string& varname) const;

  /// Return the dimensions of a variable "varname" in the specified
  /// "scope"
  virtual IntVector size(const std::string& scope, const std::string& varname) const;

  /// Return "true" if a variable exists, "false" otherwise
  virtual bool exist(const std::string& varname) const;

  /// Return "true" if the variable "varname" exists in "scope"
  virtual bool exist(const std::string& scope, const std::string& varname) const;

  /// Get a list of variable or argument names present in the file,
  /// where global scope is assumed if "scope" is omitted
  virtual bool get_names(std::list<std::string>& names,
			 const std::string& scope = std::string()) const;

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
		    int j = -1, int i = -1) const;

  /// Read a real scalar from the specified "scope" of the data file
  virtual bool read(Real& x, const std::string& scope,
		    const std::string& varname, int j = -1, int i = -1) const;
		  
  /// Read an integer scalar from data file
  virtual bool read(int& x, const std::string& varname,
		    int j = -1, int i = -1) const;

  /// Read an integer scalar from the specified "scope" of the data file
  virtual bool read(int& x, const std::string& scope,
		    const std::string& varname, int j = -1, int i = -1) const;

  /// Read a boolean from the data file
  virtual bool read(bool& x, const std::string& varname,
		    int j = -1, int i = -1) const;

  /// Read a boolean scalar from the specified "scope" of the data file
  virtual bool read(bool& x, const std::string& scope,
		    const std::string& varname, int j = -1, int i = -1) const;

  /// Read a string from the data file
  virtual bool read(std::string& s, const std::string& varname,
		    int j = -1) const;

  /// Read a string from the specified "scope" of the data file
  virtual bool read(std::string& s, const std::string& scope,
		    const std::string& varname, int j = -1) const;

  /// Read a slice of real numbers from the data file
  virtual bool read_slice(Vector& v, const std::string& varname,
		    int i = -1) const;

  /// Read a vector of real numbers from the data file
  virtual bool read(Vector& v, const std::string& varname,
		    int j = -1, int i = -1) const;

  /// Read a real vector from the specified "scope" of the data file
  virtual bool read(Vector& v, const std::string& scope,
		    const std::string& varname) const;

  /// Read a vector of integers from the data file
  virtual bool read(IntVector& v, const std::string& varname,
		    int j = -1, int i = -1) const;

  /// Read a vector of integers from the specified "scope" of the data
  /// file
  virtual bool read(IntVector& v, const std::string& scope,
		    const std::string& varname) const;

  /// Read a real matrix from the data file
  virtual bool read(Matrix& M, const std::string& varname,
		 int j = -1, int i = -1) const;

  /// Read a real matrix from the specified "scope of the data file
  virtual bool read(Matrix& M, const std::string& scope,
  		    const std::string& varname) const;

  /// Read a real 3D array from the data file
  virtual bool read(Array3D& A, const std::string& varname,
		    int j = -1, int i = -1) const;
  //  virtual bool read(Array3& A, const std::string& scope,
  //		    const std::string& varname) const;
  /// Read a real 4D array from the data file
  virtual bool read(Array4D& A, const std::string& varname,
		    int j = -1, int i = -1) const;

  /// Read the entire contents of the file to a string
  virtual bool read(std::string& s) const;

  /// Read the missing_value associated with a varname
  virtual bool read_missing_value(Real& missing_val, const std::string& varname) const;
	
  /// @}


  /// @name Miscellaneous
  /// @{

  /// Copy all the attributes of varname to a variable in the specified
  /// output file; if the last argument is omitted it will be assumed
  /// that the target variable has the same name as the source variable
  virtual bool copy_attributes(const std::string& varname,
			       OutputDataFile* output,
			       const std::string& output_varname = std::string())  
    const;

  /// Write the saved variables to an XML file
  virtual bool write_stored_variables(const std::string& filename) const {
    for (std::map<std::string,DataLogXml>::const_iterator it = log_.begin(); it != log_.end(); ++it) {
      it->second.write(filename, it->first);
    }
    return true;
  }

  /// @}

protected:
  mutable std::map<std::string, DataLogXml> log_;

  bool throw_exceptions_;
};

#endif
