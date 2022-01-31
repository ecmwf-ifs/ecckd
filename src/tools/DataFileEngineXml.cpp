/// @file      DataFileEngineXml.cpp
/// @brief     Implements the DataFileEngineXml class
/// @author    Robin J. Hogan
/// @copyright (C) Copyright 2016- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define FIXME
// Only compile this file if we have the GMV support library,
// containing the XML reading routines
#ifdef HAVE_GMVECSL
#include "XmlCAPI.h"
#include "FioCAPI.h"
#include "FioErrorCodes.h"
#include "EsaExitCodes.h"

#include "DataFileEngineXml.h"

static const int MAX_ARRAY_LENGTH = 100;

#include <string.h>
#include <cstdio>
#include <netcdf.h>

//#include <math.h>

#define CONDITIONAL_ERROR(message) \
  if (throw_exceptions_) {	   \
    ERROR << message;		   \
    THROW(XML_ERROR);			   \
  }				   \
  else {			   \
    WARNING << message;		   \
    ENDWARNING;			   \
    return false;		   \
  }

#define CHECK(operation, var) {					\
  int status = operation;					\
  if (status != XML_OK) {					\
    if (throw_exceptions_) {					\
      ERROR << "XML returned error code " << status		\
            << " when reading " << var;				\
      THROW(status);						\
    }								\
    else {							\
      return false;						\
    }								\
  }								\
}

//       WARNING << "XML returned error code " << status	
//	      << " when reading " << var;			
//      ENDWARNING;						


// When extracting an element of an array, we only request up to the
// element needed which may lead to the return value
// XML_SHORTER_OUTPUT, which is fine
#define PARTIAL_CHECK(operation, var) {				\
  int status = operation;					\
  if (status != XML_OK && status != XML_SHORTER_OUTPUT) {	\
    if (throw_exceptions_) {					\
      ERROR << "XML returned error code " << status		\
            << " when reading " << var;				\
      THROW(status);						\
    }								\
    else {							\
      WARNING << "XML returned error code " << status		\
	      << " when reading " << var;			\
      ENDWARNING;						\
      return false;						\
    }								\
  }								\
}

// Varnames coming in are separated into groups and subgroups with the
// delimiter ".", while the XML request uses the delimiter
// "/". Moreover, the XML driver requires one or two groups, so the
// following function does these conversions:
//   "a"       -> "main/a"
//   "a.b"     -> "a/b"
//   "a.b.c"   -> "a/b/c"
//   "a.b.c.d" -> "a/b/c.d" (the delimiter "." has no special meaning in the output)
static
std::string
translate_varname(const std::string& varname)
{
  std::string new_varname;
  if (varname.find_first_of(".") == std::string::npos) {
    new_varname = "main/" + varname;
  } 
  else {
    int n_nested = 0;
    new_varname = varname;
    for (size_t i = 0; i < new_varname.size(); ++i) {
      if (new_varname[i] == '.') {
	new_varname[i] = '/';
	++n_nested;
	if (n_nested > 1) {
	  break;
	}
      }
    }
  }
  //  LOG << "Name conversion: \"" << varname << "\" to \"" << new_varname << "\"\n";
  return new_varname;  
}

/// For logging configuration options we need to use a variable name
/// where "." is replaced by "_".
static
std::string
translate_varname_to_underscores(const std::string& varname)
{
  std::string new_varname = varname;
  for (size_t i = 0; i < new_varname.size(); ++i) {
    if (new_varname[i] == '.') {
      new_varname[i] = '_';
    }
  }
  return new_varname;
}


// Open the file with an absolute path
DataFileEngineXml::
DataFileEngineXml(const std::string& filename, bool throw_exceptions)
  : throw_exceptions_(throw_exceptions)
{
  //  LOG << "Reading " << filename << "\n";
  xml_load_config_file(filename.c_str());
  filename_ = filename;
}

// Close the file
DataFileEngineXml::
~DataFileEngineXml()
{ }

// Get the size of a variable
IntVector 
DataFileEngineXml::size(const std::string& varname) const
{
  IntVector out;
  FIXME
  return out;
}

IntVector 
DataFileEngineXml::size(const std::string& scope, const std::string& varname) const
{
  IntVector out;
  FIXME
  return out;    
}


// Read single real numbers
bool
DataFileEngineXml::
read(Real& x, const std::string& varname, int j, int i) const
{
  if (exist(varname)) {
  std::string new_varname = translate_varname(varname);
  if (i < 0) {
    if (j < 0) {
      double my_x;
      // initialize to all NaNs for debugging (needs include <math.h>)
      //      my_x= NAN;
      CHECK(xml_load_scalar_double(filename_.c_str(), new_varname.c_str(), &my_x),
	    new_varname);
      x = my_x;
      log_[varname].save(my_x);
    }
    else {
      int n = j+1;
      std::vector<double> my_x(n);
      // initialize to all NaNs for debugging (needs include <math.h>)
      //      for (int k=0; k<n; ++k) my_x[k] = NAN;      
      PARTIAL_CHECK(xml_load_vector_double(filename_.c_str(), new_varname.c_str(), 
				   &my_x[0], &n), new_varname);
      x = my_x[j];
      log_[varname].save(my_x);
    }
  }
  else {
    int m = i+1;
    int n = j+1;
    std::vector<double> my_x(m*n);
    // initialize to all NaNs for debugging (needs include <math.h>)
    //    for (int k=0; k<m*n; ++k) my_x[k] = NAN;      
    PARTIAL_CHECK(xml_load_matrix_double(filename_.c_str(), new_varname.c_str(), &my_x[0],
					   &n, &m), new_varname);
    x = my_x[m*j+i];
  }
  return true;
  } else {
    return false;
  }
}

bool
DataFileEngineXml::
read(Real& x, const std::string& scope, 
     const std::string& varname, int j, int i) const
{
  return read(x, scope + "." + varname, j, i);
}

// Read single integers
bool
DataFileEngineXml::
read(int& x, const std::string& varname, int j, int i) const
{
  if (exist(varname)) {
  std::string new_varname = translate_varname(varname);
  if (i < 0) {
    if (j < 0) {
      long my_x;
      CHECK(xml_load_scalar_long(filename_.c_str(), new_varname.c_str(), &my_x),
	    new_varname);
      x = my_x;
      log_[varname].save(my_x);
    }
    else {
      int n = j+1;
      std::vector<long> my_x(n);
      PARTIAL_CHECK(xml_load_vector_long(filename_.c_str(), new_varname.c_str(), 
				 &my_x[0], &n), new_varname);
      x = my_x[j];
      log_[varname].save(my_x);
    }
  }
  else {
    int m = i+1;
    int n = j+1;
    std::vector<long> my_x(m*n);
    PARTIAL_CHECK(xml_load_matrix_long(filename_.c_str(), new_varname.c_str(), &my_x[0],
			       &n, &m), new_varname);
    x = my_x[m*j+i];
  }
  return true;
  }else{
  return false;
  }
}

bool
DataFileEngineXml::
read(int& x, const std::string& scope, 
     const std::string& varname, int j, int i) const
{
  return read(x, scope + "." + varname, j, i);
}


// Read single booleans (stored as an integer in the file)
bool
DataFileEngineXml::
read(bool& x, const std::string& varname, int j, int i) const
{
  int my_x;
  bool status = read(my_x, varname, j, i);
  if (status) {
    if (my_x) {
      x = true;
    }
    else {
      x = false;
    }
  }
  return status;
}

bool
DataFileEngineXml::
read(bool& x, const std::string& scope, 
     const std::string& varname, int j, int i) const
{
  return read(x, scope + "." + varname, j, i);
}



// Read strings
bool
DataFileEngineXml::
read(std::string& s, const std::string& varname, int j) const
{
  int len = MAX_ARRAY_LENGTH-1;
  char str[MAX_ARRAY_LENGTH];
  std::string new_varname = translate_varname(varname);
  CHECK(xml_load_text(filename_.c_str(), new_varname.c_str(), str, &len),
	new_varname);
  str[len] = '\0';
  log_[varname].save(str);

  if (j > -1) {
    std::string buf = str;
    size_t istart = 0;
    for (int i = 0; i < j; ++i) {
      istart = buf.find_first_of(" \t\n", istart);
      if (istart == std::string::npos) {
	// Substring index too large
	return false;
      }
      else {
	++istart;
      }
    }
    size_t iend = buf.find_first_of(" \t\n", istart);
    if (iend != std::string::npos) {
      s = buf.substr(istart, iend-istart);
    }
    else {
      s = buf.substr(istart);
    }
  }
  else {
    s = str;
  }
  return true;
}

bool
DataFileEngineXml::
read(std::string& s, const std::string& scope,
     const std::string& varname, int j) const
{
  return read(s, scope + "." + varname, j);
}

// Read slice of real numbers - only working for NetCDF type
bool
DataFileEngineXml::
read_slice(Vector& v, const std::string& varname, int i) const
{
  CONDITIONAL_ERROR("Attempt to read vector slice for XML entry");
}


// Read vector of real numbers
bool
DataFileEngineXml::
read(Vector& v, const std::string& varname, int j, int i) const
{
  if (exist(varname)) {
  if (i >= 0) {
    CONDITIONAL_ERROR("Cannot read a vector from a 3D array in XML");
    return false;
  }
  std::string new_varname = translate_varname(varname);
  if (j < 0) {
    Vector tmp(MAX_ARRAY_LENGTH);
    int len = MAX_ARRAY_LENGTH;
    if (xml_load_vector_double(filename_.c_str(), new_varname.c_str(), 
			       tmp.data(), &len) == XML_OK) {
      v = tmp(range(0,len-1));
      log_[varname].save(v);
    }
    else {
      // We may have a vector of length one (a scalar)
      double my_x;
      CHECK(xml_load_scalar_double(filename_.c_str(), new_varname.c_str(), &my_x),
	    new_varname);
      v.resize(1);
      v[0] = my_x;
    }
  }
  else {
    int m = MAX_ARRAY_LENGTH;
    int n = j+1;
    Vector tmp(m*n);
    PARTIAL_CHECK(xml_load_matrix_double(filename_.c_str(), new_varname.c_str(), 
				 tmp.data(), &n, &m), new_varname);
    v = tmp(range(m*j,(m+1)*j-1));
  }
  return true;
    }else{
      return false;
    }
}

bool
DataFileEngineXml::
read(Vector& v, const std::string& scope,
     const std::string& varname) const
{
  return read(v, scope + "." + varname);
}

// Read vector of integers
bool
DataFileEngineXml::
read(IntVector& v, const std::string& varname, int j, int i) const
{
  if (exist(varname)){
  typedef Array<1,long> myLongVector;
  if (i >= 0) {
    CONDITIONAL_ERROR("Cannot read a vector from a 3D array in XML");
    return false;
  }
  std::string new_varname = translate_varname(varname);
  if (j < 0) {
    myLongVector tmp(MAX_ARRAY_LENGTH);
    int len = MAX_ARRAY_LENGTH;
    if (xml_load_vector_long(filename_.c_str(), new_varname.c_str(), 
			     tmp.data(), &len) == XML_OK) {
      v = tmp(range(0,len-1));
      log_[varname].save(v);
    }
    else {
      // We may have a vector of length one (a scalar)
      long my_x;
      CHECK(xml_load_scalar_long(filename_.c_str(), new_varname.c_str(), &my_x),
	    new_varname);
      v.resize(1);
      v[0] = my_x;
      log_[varname].save(v);
    }
  }
  else {
    int m = MAX_ARRAY_LENGTH;
    int n = j+1;
    myLongVector tmp(m*n);
    PARTIAL_CHECK(xml_load_matrix_long(filename_.c_str(), new_varname.c_str(), 
			       tmp.data(), &n, &m), new_varname);
    v = tmp(range(m*j,(m+1)*j-1));
  }
  return true;
  }else{
    return false;
  }
}

bool
DataFileEngineXml::
read(IntVector& v, const std::string& scope,
     const std::string& varname) const
{
  return read(v, scope + "." + varname);
}




bool
DataFileEngineXml::
read(Matrix& M, const std::string& varname, int j, int i) const
{
  if (j >= 0 || i >= 0) {
    CONDITIONAL_ERROR("Cannot read arrays with more than 2 dimensions from XML file");
    return false;
  }
  std::string new_varname = translate_varname(varname);
  int m = MAX_ARRAY_LENGTH;
  int n = MAX_ARRAY_LENGTH;
  Vector tmp(m*n);
  CHECK(xml_load_matrix_double(filename_.c_str(), new_varname.c_str(), 
			       tmp.data(), &n, &m), new_varname);
  M.resize(n,m);
  for (int k = 0; k < n; ++k) {
    M(k,__) = tmp(range(k*m,(k+1)*m-1));
  }
  return true;
}

bool
DataFileEngineXml::
read(Matrix& M, const std::string& scope,
     const std::string& varname) const
{
  return read(M, scope + "." + varname);
}



bool
DataFileEngineXml::
read(Array3D& M, const std::string& varname, int j, int i) const
{
  CONDITIONAL_ERROR("Cannot read 3D arrays from XML files");
  return false;
}


bool
DataFileEngineCfg::
read(Array4D& M, const std::string& varname, int j, int i) const
{
  CONDITIONAL_ERROR << "Cannot read 4D array from cfg file";
  return false;
}

// Return the entire contents of the file as a string
bool
DataFileEngineXml::
read(std::string& s) const
{
  DEBUG << "Cannot read the entire contents of an XML file into a string \n";
  return false;
}

// Get a list of variable or argument names present in the file. If
// "scope" is omitted then 
bool
DataFileEngineXml::
get_names(std::list<std::string>& names, const std::string& scope) const
{
  ERROR << "DataFileEngineXml::get_names not yet implemented for XML files";
  THROW(UNEXPECTED_EXCEPTION);
}


// Copy all the attributes of varname to a variable in the specified
// output file; if the last argument is omitted it will be assumed
// that the target variable has the same name as the source variable
bool
DataFileEngineXml::
copy_attributes(const std::string& varname,
		OutputDataFile* output, const std::string& output_varname) const
{
  ERROR << "DataFileEngineXml::copy_attributes not yet implemented for XML files";
  THROW(UNEXPECTED_EXCEPTION);
}


bool
DataFileEngineXml::
exist(const std::string& varname) const
{
  //CONDITIONAL_ERROR("Cannot test for existence of variables in XML files");
  std::string new_varname = translate_varname(varname);
  char a_description[MAX_ARRAY_LENGTH];
  char a_type[MAX_ARRAY_LENGTH];
  char a_units[MAX_ARRAY_LENGTH];
  int a_description_length = MAX_ARRAY_LENGTH;
  int a_type_length = MAX_ARRAY_LENGTH;
  int a_units_length = MAX_ARRAY_LENGTH;
  int a_dimensions[2];
  int status;
  status=xml_enquire(filename_.c_str(), new_varname.c_str(), 
		     a_description, &a_description_length,
		     a_type, &a_type_length,
		     a_units, &a_units_length,
		     a_dimensions);
  if (status == XML_OK || status == XML_SHORTER_OUTPUT){
    return true;
  }
  else{
    return false;
  }
}

bool
DataFileEngineXml::
exist(const std::string& scope, const std::string& varname) const
{
  return exist(scope + "." + varname);
}


bool
DataFileEngineXml::read_missing_value(Real& missing_val, const std::string& varname) const
{
  return read(missing_val, varname + "/missing_value");
}

int
DataLogXml::write(const std::string& filename, const std::string& varname) const
{
  static const char* null_str = "";
  std::string new_varname = translate_varname_to_underscores(varname);
  int status = FIO_OK;
  switch (type_) {
  case LOG_STRING:
    status = fio_write_text_sph(filename.c_str(), new_varname.c_str(), str_.c_str());
    break;
  case LOG_VECTOR:
    status = fio_write_vector_double_sph(filename.c_str(), new_varname.c_str(), &vector_[0],
					 vector_.size(), null_str, null_str);
    break;
  case LOG_INT_VECTOR:
    status = fio_write_vector_long_sph(filename.c_str(), new_varname.c_str(), &int_vector_[0],
				       int_vector_.size(), null_str, null_str);
    break;
  case LOG_REAL:
    status = fio_write_scalar_double_sph(filename.c_str(), new_varname.c_str(), real_, null_str, null_str);
    break;
  case LOG_INT:
    status = fio_write_scalar_long_sph(filename.c_str(), new_varname.c_str(), int_, null_str, null_str);
    break;
  default:
    status = FIO_UNKNOWN;
  }
  if (status != FIO_OK) {
    int len = 500;
    char my_message[500];
    fio_get_status_msg(status, my_message, &len);
    WARNING << "Unable to write config variable " << new_varname << ": " << my_message;
    ENDWARNING;
  }
  return status;
}


#endif
