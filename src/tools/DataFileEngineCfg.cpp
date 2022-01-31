/// @file      DataFileEngineCfg.cpp
/// @brief     Implements the DataFileEngineCfg class
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#include "DataFileEngineCfg.h"
#include "OutputDataFile.h"

#include <string.h>
#include <cstdio>

#define CONDITIONAL_ERROR(message) \
  if (throw_exceptions_) {	   \
    ERROR << message;		   \
    THROW(UNEXPECTED_EXCEPTION);  \
  }				   \
  else {			   \
    WARNING << message;		   \
    ENDWARNING;			   \
    return false;		   \
  }


void
DataFileEngineCfg::
open_absolute(const std::string& filename)
{
  if (data_) {
    if (!rc_append(data_, filename.c_str(), stderr)) {
      ERROR << "Error reading " << filename 
	    << " into existing config structure";
      THROW(CANNOT_OPEN_MANDATORY_FILE);	
    }
  }
  else {
    data_ = rc_read(filename.c_str(), stderr);
    if (!data_) {
      ERROR << "Error reading " << filename;
      THROW(CANNOT_OPEN_MANDATORY_FILE);
    }
  }
  DETAIL << "Reading contents of " << filename << "\n";
  // Store the filename: this way we can use DataFileEngine::name() to
  // extract it
  filename_ = filename;
}

DataFileEngineCfg::
~DataFileEngineCfg()
{
  rc_clear(data_);
  data_ = 0;
}

DataFileEngineCfg::
DataFileEngineCfg(int argc, const char** argv, bool throw_exceptions)
  : throw_exceptions_(throw_exceptions)
{
  // Extract file names from argument list, i.e. those arguments not
  // beginning with "-" or containing "=". The file names are stored
  // in the config structure with params of 0, 1 etc.
  data_ = rc_read(NULL, stderr); // Returns empty structure
  rc_register_files(data_, argc, argv);

  // Get name of first file on command line
  int ifile = rc_get_file(argc, argv);

  if (ifile) {
    // Read configuration information from the file.
    open_absolute(argv[ifile]);
  }
  // Supplement configuration information with command-line arguments
  rc_register_args(data_, argc, argv);    
}

// Get the size of a variable
IntVector 
DataFileEngineCfg::size(const std::string& varname) const
{
  IntVector out;
  int m, n, len;
  len = rc_size(data_, varname.c_str(), &m, &n);
  if (m == 0 && n == 0) {
    if (len > 0) {
      out.resize(1);
      out = len;
    } // Else out will be empty indicating varname does not exist or
      // is empty
  }
  else {
    if (m > 0 && n > 0) {
      out.resize(2);
      out(0) = m;
      out(1) = n;
      if (len != m*n) {
	WARNING << "Matrix \"" << varname << "\" in cfg file \"" << filename_
		<< "\" has length inconsistent with dimensions ("
		<< len << " != " << m << " * " << n << ")";
	WARNING << rc_sprint(data_);
	ENDWARNING;
      }
    }
    else {
      out.resize(1);
      out = len;
      if (len != m+n) {
	WARNING << "Stated length of vector in cfg file is not correct ("
		<< len << " != " << m+n << ")";
	ENDWARNING;
      }
    }
  }
  return out; 
}

IntVector 
DataFileEngineCfg::size(const std::string& scope, const std::string& varname) const
{
  IntVector out;
  rc_set_section(data_, scope.c_str());
  out = size(varname);
  rc_set_section(data_, 0);
  return out;
}


// Read single real numbers
bool
DataFileEngineCfg::
read(Real& x, const std::string& varname, int j, int i) const
{
  if (j >= 0 || i >= 0) {
    if (j < 0) j = 0;
    if (i < 0) i = 0;
    if (j <= 0 || i <= 0) {
      return rc_assign_real_element(data_, varname.c_str(), i+j, &x);
    }
    else {
      int m, n;
      rc_real** data = rc_get_real_matrix(data_, varname.c_str(), &m, &n);
      if (!data) {
	return false;
      }
      if (j > m-1 || i < n-1) {
	CONDITIONAL_ERROR("Indices to cfg matrix out of range");
	rc_free(data);
	return false;
      }
      x = data[j][i];
      rc_free(data);
      return true;
    }
  }
  else {
    return rc_assign_real(data_, varname.c_str(), &x);
  }
  return false;
}

bool
DataFileEngineCfg::
read(Real& x, const std::string& scope, 
     const std::string& varname, int j, int i) const
{
  rc_set_section(data_, scope.c_str());
  bool status = read(x, varname, j, i);
  rc_set_section(data_, 0);
  return status;
}

// Read single integers
bool
DataFileEngineCfg::
read(int& x, const std::string& varname, int j, int i) const
{
  if (j >= 0 || i >= 0) {
    CONDITIONAL_ERROR("Attempt to load specific integer from cfg file");
  }
  return rc_assign_int(data_, varname.c_str(), &x);
}

bool
DataFileEngineCfg::
read(int& x, const std::string& scope, 
     const std::string& varname, int j, int i) const
{
  rc_set_section(data_, scope.c_str());
  bool status = rc_assign_int(data_, varname.c_str(), &x);
  rc_set_section(data_, 0);
  return status;
}


// Read single booleans
bool
DataFileEngineCfg::
read(bool& x, const std::string& varname, int j, int i) const
{
  if (j >= 0 || i >= 0) {
    CONDITIONAL_ERROR("Attempt to load specific boolean from cfg file");
  }
  x = rc_get_boolean(data_, varname.c_str());
  return true;
}

bool
DataFileEngineCfg::
read(bool& x, const std::string& scope, 
     const std::string& varname, int j, int i) const
{
  rc_set_section(data_, scope.c_str());
  x = rc_get_boolean(data_, varname.c_str());
  rc_set_section(data_, 0);
  return true;

}


// Read strings
bool
DataFileEngineCfg::
read(std::string& s, const std::string& varname, int j) const
{
  if (j >= 0) {
    char* str = rc_get_substring(data_, varname.c_str(), j);
    if (!str) {
      return false;
    }
    else {
      s = str;
      rc_free(str);
      return true;
    }
  }
  else {
    char* str = rc_get_string(data_, varname.c_str());
    if (!str) {
      return false;
    }
    else {
      s = str;
      rc_free(str);
      return true;
    }
  }
  return false;
}

bool
DataFileEngineCfg::
read(std::string& s, const std::string& scope,
     const std::string& varname, int j) const
{
  rc_set_section(data_, scope.c_str());
  bool status = read(s, varname, j);
  rc_set_section(data_, 0);
  return status;
}

// Read slice of real numbers - only working for NetCDF type
bool
DataFileEngineCfg::
read_slice(Vector& v, const std::string& varname, int i) const
{
  CONDITIONAL_ERROR("Attempt to read vector slice for cfg entry");
}

// Read vector of real numbers
bool
DataFileEngineCfg::
read(Vector& v, const std::string& varname, int j, int i) const
{
  bool is_matrix = rc_is_matrix(data_, varname.c_str());

  if (i >= 0) {
    CONDITIONAL_ERROR("Attempt to read cfg entry as if it was a 3D array");
  }
  else if (j > 0 && !is_matrix) {
    CONDITIONAL_ERROR("Attempt to read a row from a cfg entry that is not a matrix");
  }
  else if (is_matrix) {
    int m, n;
    rc_real** data = rc_get_real_matrix(data_, varname.c_str(), &m, &n);
    if (!data) {
      return false;
    }
    
    Matrix M_tmp(*data, expression_size(m, n));
    v.resize(n);
    v = M_tmp(j,__);
    
    rc_free_matrix(reinterpret_cast<void**>(data));
    return true;
    
  }
  else {
    int len;
    rc_real* data = rc_get_real_vector(data_, varname.c_str(), &len);
    if (!data) {
      return false;
    }
    v.resize(len);
    Vector v_tmp(data,expression_size(len));
    v = v_tmp;
    rc_free(data);
    return true;
  }
  return false;
}

bool
DataFileEngineCfg::
read(Vector& v, const std::string& scope,
     const std::string& varname) const
{
  rc_set_section(data_, scope.c_str());
  bool status = read(v, varname);
  rc_set_section(data_, 0);
  return status;    
}

// Read vector of integers
bool
DataFileEngineCfg::
read(IntVector& v, const std::string& varname, int j, int i) const
{
  bool is_matrix = rc_is_matrix(data_, varname.c_str());

  if (i >= 0) {
    CONDITIONAL_ERROR("Attempt to read cfg entry as if it was a 3D array");
  }
  else if (j > 0 && !is_matrix) {
    CONDITIONAL_ERROR("Attempt to read a row from a cfg entry that is not a matrix");
  }
  else if (is_matrix) {
    int m, n;
    int** data = rc_get_int_matrix(data_, varname.c_str(), &m, &n);
    if (!data) {
      return false;
    }
    
    IntMatrix M_tmp(*data,expression_size(m, n));
    v.resize(n);
    v = M_tmp(j,__);
    
    rc_free_matrix(reinterpret_cast<void**>(data));
    //      rc_free_matrix(data);
    return true;
    
  }
  else {
    int len;
    int* data = rc_get_int_vector(data_, varname.c_str(), &len);
    if (!data) {
      return false;
    }
    v.resize(len);
    IntVector v_tmp(data,expression_size(len));
    v = v_tmp;
    rc_free(data);
    return true;
  }
  
  return false;
}

bool
DataFileEngineCfg::
read(IntVector& v, const std::string& scope,
     const std::string& varname) const
{
  rc_set_section(data_, scope.c_str());
  bool status = read(v, varname);
  rc_set_section(data_, 0);
  return status;    
}




bool
DataFileEngineCfg::
read(Matrix& M, const std::string& varname, int j, int i) const
{
  if (j >= 0 || i >= 0) {
    CONDITIONAL_ERROR("Attempt to read sub-part of cfg matrix");
  }
  int m, n;
  rc_real** data = rc_get_real_matrix(data_, varname.c_str(), &m, &n);
  if (!data) {
    return false;
  }
  M.resize(m, n);
  Matrix M_tmp(*data,expression_size(m, n));
  M = M_tmp;
  rc_free(data);
  return true;
}

bool
DataFileEngineCfg::
read(Matrix& M, const std::string& scope,
     const std::string& varname) const
{
  rc_set_section(data_, scope.c_str());
  bool status = read(M, varname);
  rc_set_section(data_, 0);
  return status;    
}



bool
DataFileEngineCfg::
read(Array3D& M, const std::string& varname, int j, int i) const
{
    ERROR << "Cannot read 3D array from cfg file";
    THROW(READ_ERROR_PRODUCT_MODEL);
}

bool
DataFileEngineCfg::
read(Array4D& M, const std::string& varname, int j, int i) const
{
    ERROR << "Cannot read 4D array from cfg file";
    THROW(READ_ERROR_PRODUCT_MODEL);
}


// Return the entire contents of the file as a string
bool
DataFileEngineCfg::
read(std::string& s) const
{
  char* str = rc_sprint(data_);
  if (str) {
    s = str;
    rc_free(str);
    return true;
  }
  else {
    return false;
  }
}

// Get a list of variable or argument names present in the file. If
// "scope" is omitted then 
bool
DataFileEngineCfg::
get_names(std::list<std::string>& names, const std::string& scope) const
{
  if (scope.empty()) {
    // Return all variables
    rc_data* data = data_;
    names.clear();
    while (data && data->param) {
      names.push_back(data->param);
      data = data->next;
    }
  }
  else if (scope == DATA_FILE_GLOBAL_SCOPE) {
    WARNING << "CFG files do not have the concept of global attributes";
    ENDWARNING;
    return 0;
  }
  else {
    // Return variables prefixed by contents of "scope"
    rc_data* data = data_;
    std::string param;
    std::string scope_dot = scope + ".";
    size_t scope_len = scope_dot.size();
    names.clear();
    while (data && data->param) {
      param = data->param;
      if (param.size() > scope_len
	  && param.substr(0,scope_len) == scope_dot) {
	names.push_back(param.substr(scope_len, std::string::npos));
      }
      data = data->next;
    }
  }
  
  return true;
}


// Copy all the attributes of varname to a variable in the specified
// output file; if the last argument is omitted it will be assumed
// that the target variable has the same name as the source variable
bool
DataFileEngineCfg::
copy_attributes(const std::string& varname,
		OutputDataFile* output, const std::string& output_varname) const
{
  std::string out_varname;
  if (output_varname.empty()) {
    out_varname = varname;
  }
  else {
    out_varname = output_varname;
  }
  if (!output->exist(out_varname)) {
    return false;
  }

  rc_data* data = data_;
  while ((data = rc_find_next_simple_member(data, varname.c_str()))) {
    const char* member_name = data->param+varname.size()+1;
    // First try to read this as a vector of numbers
    int len;
    rc_set_section(data, varname.c_str());
    rc_real* x = rc_get_real_vector(data, member_name, &len);
    rc_pop_section(data);
    if (!x) {
      // Assume this is a string
      output->write(data->value, out_varname, member_name);
    }
    else {
      Vector v(x, expression_size(len));
      rc_free(x);
      output->write(v, out_varname, member_name);
    }
    data = data->next;
  }
  
  return true;
}


bool
DataFileEngineCfg::
exist(const std::string& varname) const
{
  return rc_exists(data_, varname.c_str());
}

bool
DataFileEngineCfg::
exist(const std::string& scope, const std::string& varname) const
{
  rc_set_section(data_, scope.c_str());
  bool status = rc_exists(data_, varname.c_str());
  rc_set_section(data_, 0);
  return status;
}


bool
DataFileEngineCfg::read_missing_value(Real& missing_val, const std::string& varname) const
{
  rc_set_section(data_, varname.c_str());
  bool status = rc_assign_real(data_, "missing_value", &missing_val);
  rc_set_section(data_, 0);
  return status;
}
