/// @file      OutputDataFile.cpp
/// @brief     Implements the OutputDataFile class
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#include <cstdio>
#include <ctime>
#include <netcdf.h>

//#include "State.h"
#include "OutputDataFile.h"

std::string OUTPUT_TIME_DIMENSION;
std::string OUTPUT_RANGE_DIMENSION;
static int ncstatus;
#define NC_CHECK_UNNAMED(a) if ((ncstatus = (a)) != NC_NOERR) {	\
    ERROR << nc_strerror(ncstatus) << " [NC code " << ncstatus << "]"; \
    THROW(UNEXPECTED_EXCEPTION); }
#define NC_CHECK(a, name) if ((ncstatus = (a)) != NC_NOERR) {	\
    ERROR << nc_strerror(ncstatus) << " [NC code " << ncstatus << "]" \
	  << ": " << name; THROW(UNEXPECTED_EXCEPTION); }
//#define NC_CHECK(a, name) if ((ncstatus = (a)) != NC_NOERR) {		
//    WARNING << nc_strerror(ncstatus) << ": " << name; ENDWARNING; }

#define CONDITIONAL_ERROR(message) \
  if (throw_exceptions_) {	   \
    ERROR << message;		   \
    THROW(UNEXPECTED_EXCEPTION);			   \
  }				   \
  else {			   \
    WARNING << message;		   \
    ENDWARNING;			   \
    return false;		   \
  }

static
nc_type
type2nctype(OutputDataFileVariableType type) {
  if (type == DOUBLE) {
    return NC_DOUBLE;
  }
  else if (type == FLOAT) {
    return NC_FLOAT;
  }
  else if (type == INT) {
    return NC_INT;
  }
  else if (type == SHORT) {
    return NC_SHORT;
  }
  else if (type == STRING) {
    return NC_CHAR;
  }
  else if (type == BYTE) {
    return NC_BYTE;
  }
  else {
    ERROR << "Type not recognized";
    THROW(PRODUCT_FORMAT_ERROR);
  }
}


void
OutputDataFile::
open_absolute(const std::string& filename, OutputDataFileMode mode, 
	      DataFileType type)
{
  /* Close existing file first */
  if (type_ != UNDETERMINED) {
    close();
  }
  is_netcdf4hdf5_ = false;
  if (type == UNDETERMINED) {
    // Determine type from filename
    size_t dotpos = filename.find_last_of('.');
    if (dotpos == std::string::npos) {
      ERROR << "Cannot determine file type for \"" << filename << "\"";
      THROW(PRODUCT_FORMAT_ERROR);
    }
    std::string extension = filename.substr(dotpos+1);
    if (extension == "nc" || extension == "NC"
	|| extension == "cdf" || extension == "CDF") {
      type = NETCDF;
      //#ifdef NC_NETCDF4
      //      is_netcdf4hdf5_ = true;
      //#endif
    }
    else if (extension == "h5" || extension == "hdf") {
#ifndef NC_NETCDF4
      ERROR << "Cannot write NetCDF-4/HDF-5 file \"" << filename 
	    << "\": compiled with NetCDF-3 library";
      THROW(PRODUCT_FORMAT_ERROR);
#else
      type = NETCDF;
      is_netcdf4hdf5_ = true;
#endif
    }
    else if (extension == "cfg") {
      type = CONFIG;
    }
    else {
      ERROR << "File extension \"" << extension << "\" not recognised";
      THROW(PRODUCT_FORMAT_ERROR);
    }
  }


  if (type == NETCDF) {
    // nc_mode is a bit-field of options
    int nc_mode = 0;
    
#ifdef NC_NETCDF4
    if (is_netcdf4hdf5_) {
      nc_mode |= NC_NETCDF4;
    }
#endif
    if (mode == OUTPUT_MODE_NOCLOBBER || mode == OUTPUT_MODE_APPEND) {
      nc_mode |= NC_NOCLOBBER;
    }
    /*
    if (is_netcdf4hdf5_) {
      // Check if we have NetCDF-4 capability
#ifdef NC_NETCDF4
      if (mode == OUTPUT_MODE_APPEND) {
	DETAIL << "Opening " << filename << " for appending (NetCDF-4 format)\n";
	NC_CHECK(nc_open(filename.c_str(), NC_WRITE, &root_ncid_), filename);
	NC_CHECK(nc_inq_ncid(root_ncid_, "HeaderData", &header_ncid_), filename);
	NC_CHECK(nc_inq_ncid(root_ncid_, "ScienceData", &ncid_), filename);
      }
      else {
	DETAIL << "Opening " << filename << " for writing (NetCDF-4 format)\n";
	// split data into 2 groups, ScienceData and Header
	NC_CHECK(nc_create(filename.c_str(), nc_mode, &root_ncid_), filename);
	NC_CHECK(nc_def_grp(root_ncid_, "HeaderData", &header_ncid_), filename);
	NC_CHECK(nc_def_grp(root_ncid_, "ScienceData", &ncid_), filename);
	LOG << "  Creating two NetCDF-4 groups: HeaderData and ScienceData\n";
      }
#else
      ERROR << "Compiled with NetCDF3: unable to create groups";
      THROW(WRITE_ERROR);
#endif
    }
    else */ {
      if (mode == OUTPUT_MODE_APPEND) {
	ERROR << "Cannot append a NetCDF-3 file";
	THROW(WRITE_ERROR);
      }
      DETAIL << "Opening " << filename << " for writing (NetCDF-3 format)\n";
      NC_CHECK(nc_create(filename.c_str(), nc_mode, &root_ncid_), filename);
      ncid_=root_ncid_;
    }
  }
  else if (type == CONFIG) {
    DETAIL << "Opening " << filename << " for writing\n";
    file_ = fopen(filename.c_str(), "w");
    if (!file_) {
      ERROR << "Error writing " << filename;
      THROW(WRITE_ERROR);
    }
  }
  type_ = type;
  filename_ = filename;
}

void
OutputDataFile::
open(const std::string& filename, OutputDataFileMode mode, DataFileType type)
{
  if (filename[0] == '/' || filename[0] == '.') {
    open_absolute(filename, mode, type);
  }
  else {
    std::string absfilename;
    absfilename = first_directory();
    if (!absfilename.empty()) {
      absfilename += "/";      
    }
    absfilename += filename;
    open_absolute(absfilename, mode, type);
  }
}

void
OutputDataFile::
close()
{
  if (type_ == CONFIG) {
    fclose(file_);
    file_ = 0;
    LOG << "Closed " << filename_ << "\n";
  }
  else if (type_ == NETCDF) {
    NC_CHECK(nc_close(root_ncid_), filename_);
    root_ncid_ = 0;
    LOG << "Closed " << filename_ << "\n";
  }
  type_ = UNDETERMINED;
}

void
OutputDataFile::
set_options(DataFile& config)
{
  if (type_ == NETCDF) {
    //set here the default for the output_time_dimension name. 
    //default is "time" and it can be externally changed
    //see there the output_height_dimension name. In HDF5
    //the dimension cannot be a variable
    if (config.read(time_dimension_name_, "time_dimension_name")){
      OUTPUT_TIME_DIMENSION = time_dimension_name_;
      
    }
    else{
      OUTPUT_TIME_DIMENSION = "time";
    }
    if (config.read(range_dimension_name_, "range_dimension_name")){
      OUTPUT_RANGE_DIMENSION = range_dimension_name_; 
    }
    else{
      OUTPUT_RANGE_DIMENSION = "height";
    }
    LOG << "output file horizontal dimension will be " << OUTPUT_TIME_DIMENSION << "\n";
    LOG << "output file vertical dimension will be " << OUTPUT_RANGE_DIMENSION << "\n";

    // Option for reducing the size of output file by removing the
    // avg kernel and error_correlation variables
    // defalut: full error statistics
    //write_full_err_stats_ = 1;
    if (!config.read(write_full_err_stats_, "write_full_err_stats") || write_full_err_stats_ == 1){
      LOG << "all error statistics will be written to output file \n";
    }else{
      LOG << "reduced error statistics in output file \n";
    }

    // Options for chunking and compression of NetCDF-4/HDF-5 files
    if (config.read(chunk_size_, "chunk_size") && is_netcdf4hdf5_){
      if (chunk_size_ == -1) {
	LOG << "  NetCDF has horizontal dimension unlimited (default)\n";
      }
      else if (chunk_size_ == 0) {
	LOG << "  NetCDF has plain dimensions\n";
      }
      else if (chunk_size_ > 0) {
	LOG << "  NetCDF file will be organized into " << chunk_size_ << " chunks for the horizontal dimension\n";
      }
    }
    else {
      LOG << "  No chunking performed on NetCDF. Time unlimited instead.\n";
    }
  }
}


void
OutputDataFile::
define_dimension(const std::string& dimname, int length)
{
  // This function does nothing if the type is not netcdf
  if (type_ == NETCDF) {
    int dimid;
    if (length > 0) {
      int dimid;
      // Only create dimension if it doesn't yet exist
      if (nc_inq_dimid(ncid_, dimname.c_str(), &dimid) != NC_NOERR) {
	NC_CHECK(nc_def_dim(ncid_, dimname.c_str(), length, &dimid), dimname);
      }
    }
    else {
      NC_CHECK(nc_def_dim(ncid_, dimname.c_str(), NC_UNLIMITED, &dimid), dimname);
    }
  }
}

void
OutputDataFile::
define_variable(const std::string& varname, 
		OutputDataFileVariableType type,
		const std::string& dim1name,
		const std::string& dim2name,
		const std::string& dim3name,
		const std::string& dim4name) 
{
  // This function does nothing if the type is not netcdf
  if (type_ == NETCDF) {
    int ndims = 0;
    int dimid[4];
    if (!dim1name.empty()) {
      ndims = 1;
      NC_CHECK(nc_inq_dimid(ncid_, dim1name.c_str(), &dimid[0]), dim1name);
      if (!dim2name.empty()) {
	ndims = 2;
	NC_CHECK(nc_inq_dimid(ncid_, dim2name.c_str(), &dimid[1]), dim2name);
	if (!dim3name.empty()) {
	  ndims = 3;
	  NC_CHECK(nc_inq_dimid(ncid_, dim3name.c_str(), &dimid[2]), dim3name);
	  if (!dim4name.empty()) {
	    ndims = 4;
	    NC_CHECK(nc_inq_dimid(ncid_, dim4name.c_str(), &dimid[3]), dim4name);
	  }
	}
      }
    }
    int varid;
    NC_CHECK(nc_def_var(ncid_, varname.c_str(), type2nctype(type),
			ndims, dimid, &varid), varname);
#ifdef NC_NETCDF4
    // Chunking only the first dimension (time)
    if (chunk_size_ > 0 && is_netcdf4hdf5_) {
      int dimid;
      size_t chunks[ndims];
      std::fill(chunks,chunks+ndims,1);
      if (ndims >= 1 && dim1name == OUTPUT_TIME_DIMENSION) {
	// Check length of this dimension
	NC_CHECK(nc_inq_dimid (ncid_, dim1name.c_str(), &dimid), dim1name);
	NC_CHECK(nc_inq_dimlen(ncid_, dimid, &chunks[0]), dim1name);
	if (chunks[0] > chunk_size_) {
	  // Length is greater than chunk size so we can do chunking
	  chunks[0] = chunk_size_;
	  //DETAIL << "  Organizing the 1st dimension of " << varname << " in " << chunk_size_ << " chunks\n";
	  // Set the remaining chunk sizes to the dimension lengths
	  if (ndims >= 2) {
	    NC_CHECK(nc_inq_dimid (ncid_, dim2name.c_str(), &dimid), dim2name);
	    NC_CHECK(nc_inq_dimlen(ncid_, dimid, &chunks[1]), dim2name);
	    if (ndims >= 3) {
	      NC_CHECK(nc_inq_dimid (ncid_, dim3name.c_str(), &dimid), dim3name);
	      NC_CHECK(nc_inq_dimlen(ncid_, dimid, &chunks[2]), dim3name);
	      if (ndims >= 4) {
		NC_CHECK(nc_inq_dimid (ncid_, dim4name.c_str(), &dimid), dim4name);
		NC_CHECK(nc_inq_dimlen(ncid_, dimid, &chunks[3]), dim4name);
	      }
	    }
	  }
	  NC_CHECK(nc_def_var_chunking(ncid_, varid, NC_CHUNKED, chunks), varname);
	}
      }
    }
#endif
  }
}

void
OutputDataFile::
deflate_variable(const std::string& varname)
{
  // This function does nothing if the type is not netcdf4
  if (type_ == NETCDF) {
#ifdef NC_NETCDF4
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    // Turn on level-2 compression as well as byte shuffle
    NC_CHECK(nc_def_var_deflate(ncid_, varid, 1, 1, 2), varname);
#endif
  }
}

void
OutputDataFile::
set_chunking(const IntVector& chunk_size, const std::string& varname)
{
  // This function does nothing if the type is not netcdf4
  if (type_ == NETCDF) {
#ifdef NC_NETCDF4
    int varid;
    size_t chunks[chunk_size.size()];
    for (int ii = 0; ii < chunk_size.size(); ++ii) {
      chunks[ii] = chunk_size(ii);
    }
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    NC_CHECK(nc_def_var_chunking(ncid_, varid, NC_CHUNKED, chunks), varname);
#endif
  }
}

bool
OutputDataFile::
write_missing_value(Real missing_value,
		    const std::string& varname)
{
  if (type_ == CONFIG) {
    std::string missing_value_name = varname + ".missing_value";
    return rc_write_reals(file_, missing_value_name.c_str(), 1,
			  &missing_value);
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    nc_type type;
    NC_CHECK(nc_inq_vartype(ncid_, varid, &type), varname);

#ifdef REAL_IS_FLOAT
    NC_CHECK(nc_put_att_float(ncid_, varid, "_FillValue", 
			      type, 1, &missing_value), varname);
    NC_CHECK(nc_put_att_float(ncid_, varid, "missing_value", 
			      type, 1, &missing_value), varname);
#else
    NC_CHECK(nc_put_att_double(ncid_, varid, "_FillValue", 
			       type, 1, &missing_value), varname);
    NC_CHECK(nc_put_att_double(ncid_, varid, "missing_value", 
			       type, 1, &missing_value), varname);
#endif
    return true;
  }
  return false;
}

void
OutputDataFile::
end_define_mode()
{
  // This function does nothing if the type is not netcdf
  if (type_ == NETCDF) {
    NC_CHECK(nc_enddef(ncid_), filename_);
  }
}

bool
OutputDataFile::
write_comment(const std::string& comment)
{
  // This function does nothing if the type is not cfg
  if (type_ == CONFIG) {
    return rc_write_comment(file_, comment.c_str());
  }
  return false;
}



// Write single real numbers
bool
OutputDataFile::
write(Real x, const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write specific real number to cfg file");
    }
    return rc_write_reals(file_, varname.c_str(), 1, &x);
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid_, varid, &ndims), varname);
    if (ndims > 2) {
      CONDITIONAL_ERROR("Cannot write single element to netcdf file with more than 2 dimensions");
    }
    if (j == -1) j = 0;
    if (i == -1) i = 0;
    size_t index[2] = {static_cast<size_t>(j), static_cast<size_t>(i) };
#ifdef REAL_IS_FLOAT
    NC_CHECK(nc_put_var1_float(ncid_, varid, index, &x), varname);
#else
    NC_CHECK(nc_put_var1_double(ncid_, varid, index, &x), varname);
#endif
    return true;
  }
  return false;
}

bool
OutputDataFile::
write(Real x, const std::string& scope, 
     const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    std::string fullvarname = scope + "." + varname;
    return rc_write_reals(file_, fullvarname.c_str(), 1, &x);
  }
  else if (type_ == NETCDF) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write specific real element of a netcdf attribute");
      return false;
    }
    int varid;
    if (scope == DATA_FILE_GLOBAL_SCOPE) {
      varid = NC_GLOBAL;
    }
    else {
      NC_CHECK(nc_inq_varid(ncid_, scope.c_str(), &varid), scope);
    }
    nc_type type;
    NC_CHECK(nc_inq_vartype(ncid_, varid, &type), varname);
    // If the variable type is "double", then write attributes with
    // the same type, otherwise use "float"
    if (type != NC_DOUBLE) {
      type = NC_FLOAT;
    }
#ifdef REAL_IS_FLOAT
    NC_CHECK(nc_put_att_float(ncid_, varid, varname.c_str(), 
			      type, 1, &x), varname);
#else
    NC_CHECK(nc_put_att_double(ncid_, varid, varname.c_str(), 
			       type, 1, &x), varname);
#endif
    return true;
  }
  return false;
}

// Write single integers
bool
OutputDataFile::
write(int x, const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to load specific integer from cfg file");
    }
    return rc_write_ints(file_, varname.c_str(), 1, &x);
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid_, varid, &ndims), varname);
    if (ndims > 2) {
      CONDITIONAL_ERROR("Cannot write single element to netcdf file with more than 2 dimensions");
    }
    if (j == -1) j = 0;
    if (i == -1) i = 0;
    size_t index[2] = {static_cast<size_t>(j), static_cast<size_t>(i) };
    NC_CHECK(nc_put_var1_int(ncid_, varid, index, &x), varname);
    return true;
  }
  return false;
}

bool
OutputDataFile::
write(int x, const std::string& scope, 
     const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    std::string fullvarname = scope + "." + varname;
    return rc_write_ints(file_, fullvarname.c_str(), 1, &x);
  }
  else if (type_ == NETCDF) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write specific integer of a netcdf attribute");
      return false;
    }
    int varid;
    if (scope == DATA_FILE_GLOBAL_SCOPE) {
      varid = NC_GLOBAL;
    }
    else {
      NC_CHECK(nc_inq_varid(ncid_, scope.c_str(), &varid), scope);
    }
    NC_CHECK(nc_put_att_int(ncid_, varid, varname.c_str(), NC_INT, 1, &x),
	     varname);
    return true;
  }
  return false;
}



// Write strings
bool
OutputDataFile::
write(const std::string& s, const std::string& varname, int j)
{
  if (type_ == CONFIG) {
    if (j >= 0) {
      CONDITIONAL_ERROR("Attempt to write substring to a cfg file");
    }
    else {
      return rc_write_string(file_, varname.c_str(), s.c_str());
    }
  }
  else if (type_ == NETCDF) {
    if (j >= 0) {
      CONDITIONAL_ERROR("Attempt to write substring to a netcdf file");
    }
    else { 
      // CONDITIONAL_ERROR("Writing of netcdf string variables not yet implemented");
      NC_CHECK(nc_put_att_text(root_ncid_, NC_GLOBAL, varname.c_str(),
			       s.size(), s.c_str()), varname);
    }
  }
  return false;
}

bool
OutputDataFile::
write(const std::string& s, const std::string& scope,
     const std::string& varname, int j)
{
  if (type_ == CONFIG) {
    if (j >= 0) {
      CONDITIONAL_ERROR("Attempt to write substring to a cfg file");
    }
    std::string fullvarname = scope + "." + varname;
    return rc_write_string(file_, fullvarname.c_str(), s.c_str());
  }
  else if (type_ == NETCDF) {
    if (j >= 0) {
      CONDITIONAL_ERROR("Attempt to write substring of netcdf attribute");
    }
    int varid;
    int ncid_local;
    if (scope == DATA_FILE_GLOBAL_SCOPE) {
      varid = NC_GLOBAL;
      ncid_local=root_ncid_;
    }
    else {
      ncid_local=ncid_;
      NC_CHECK(nc_inq_varid(ncid_local, scope.c_str(), &varid), scope);
    }
    size_t len = s.size();
    NC_CHECK(nc_put_att_text(ncid_local, varid, varname.c_str(),
			     len, s.c_str()), varname);
    return true;
  }
  return false;
}

// Write vector of real numbers
bool
OutputDataFile::
write(const Vector& v, const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write sub-part of cfg vector");
    }
    return rc_write_reals(file_, varname.c_str(), v.size(),
			  v.data_pointer());
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid_, varid, &ndims), varname);
    if (ndims > 3) {
      CONDITIONAL_ERROR("Cannot write to netcdf variable with more than 3 dimensions");
    }
    else if (i < 0 && ndims == 3) {
      CONDITIONAL_ERROR("No second index specified for vector written to 3-D netcdf array");
    }
    else if (j < 0 && ndims > 1) {
      CONDITIONAL_ERROR("No index specified for vector written to 2-D or 3-D netcdf array");
    }

    int dimid[3];
    NC_CHECK(nc_inq_vardimid(ncid_, varid, dimid), varname);
    size_t len[3];
    for (int idim = 0; idim < ndims; idim++) {
      NC_CHECK_UNNAMED(nc_inq_dimlen(ncid_, dimid[idim], len+idim));
    }
    size_t start[3] = {0, 0, 0};
    if (j >= 0) { 
      start[0] = j;
      len[0] = 1;
      if (i >= 0) {
	start[1] = i;
	len[1] = 1;
      }
    }
    /*
    size_t length = 1;
    if (ndims > 0) {
      length = len[ndims-1];
    }
    */

#ifdef REAL_IS_FLOAT
    NC_CHECK(nc_put_vara_float(ncid_, varid, start, len, v.data_pointer()),
	     varname);
#else
    NC_CHECK(nc_put_vara_double(ncid_, varid, start, len, v.data_pointer()),
	     varname);
#endif
    return true;
  }
  return false;
}

bool
OutputDataFile::
write(const Vector& v, const std::string& scope,
      const std::string& varname)
{
  if (type_ == CONFIG) {
    std::string fullvarname = scope + "." + varname;
    return rc_write_reals(file_, fullvarname.c_str(),
			  v.size(), v.data_pointer());
  }
  else if (type_ == NETCDF) {
    int varid;
    if (scope == DATA_FILE_GLOBAL_SCOPE) {
      varid = NC_GLOBAL;
    }
    else {
      NC_CHECK(nc_inq_varid(ncid_, scope.c_str(), &varid), scope);
    }
    nc_type type;
    NC_CHECK(nc_inq_vartype(ncid_, varid, &type), varname);
    if (type != NC_DOUBLE) {
      type = NC_FLOAT;
    }

#ifdef REAL_IS_FLOAT
    NC_CHECK(nc_put_att_float(ncid_, varid, varname.c_str(), 
			      type, v.size(), v.data_pointer()),
	     varname);
#else
    NC_CHECK(nc_put_att_double(ncid_, varid, varname.c_str(), 
			       type, v.size(), v.data_pointer()),
	     varname);
#endif
    return true;
  }
  return false;
}
 
// Write vector of integers
bool
OutputDataFile::
write(const IntVector& v, const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write sub-part of cfg vector");
    }
    return rc_write_ints(file_, varname.c_str(), v.size(),
			 v.data_pointer());
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid_, varid, &ndims), varname);
    if (ndims > 3) {
      CONDITIONAL_ERROR("Cannot write to netcdf variable with more than 3 dimensions");
    }
    else if (i < 0 && ndims == 3) {
      CONDITIONAL_ERROR("No second index specified for vector written to 3-D netcdf array");
    }
    else if (j < 0 && ndims > 1) {
      CONDITIONAL_ERROR("No index specified for vector written to 2-D or 3-D netcdf array");
    }

    int dimid[3];
    NC_CHECK(nc_inq_vardimid(ncid_, varid, dimid), varname);
    size_t len[3];
    for (int idim = 0; idim < ndims; idim++) {
      NC_CHECK_UNNAMED(nc_inq_dimlen(ncid_, dimid[idim], len+idim));
    }
    size_t start[3] = {0, 0, 0};
    if (j >= 0) { 
      start[0] = j;
      len[0] = 1;
      if (i >= 0) {
	start[1] = i;
	len[1] = 1;
      }
    }
    /*
    size_t length = 1;
    if (ndims > 0) {
      length = len[ndims-1];
    }
    */

    NC_CHECK(nc_put_vara_int(ncid_, varid, start, len, v.data_pointer()),
	     varname);
    return true;
  }
  return false;
}
  

bool
OutputDataFile::
write(const IntVector& v, const std::string& scope,
     const std::string& varname)
{
  if (type_ == CONFIG) {
    std::string fullvarname = scope + "." + varname;
    return rc_write_ints(file_, fullvarname.c_str(),
			 v.size(), v.data_pointer());
  }
  else if (type_ == NETCDF) {
    int varid;
    if (scope == DATA_FILE_GLOBAL_SCOPE) {
      varid = NC_GLOBAL;
    }
    else {
      NC_CHECK(nc_inq_varid(ncid_, scope.c_str(), &varid), scope);
    }

    NC_CHECK(nc_put_att_int(ncid_, varid, varname.c_str(), 
			    NC_INT, v.size(), v.data_pointer()),
	     varname);
    return true;
  }
  return false;
}
  

bool
OutputDataFile::
write(const Matrix& M, const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write sub-part of cfg matrix");
    }
    CONDITIONAL_ERROR("Writing matrices to cfg files is not yet implemented");
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid_, varid, &ndims), varname);
    if (ndims > 4) {
      CONDITIONAL_ERROR("Cannot write to netcdf array with more than 4 dimensions");
    }
    else if (i < 0 && ndims == 4) {
      CONDITIONAL_ERROR("No second index specified for matrix written to 4-D netcdf array");
    }
    else if (j < 0 && ndims > 2) {
      CONDITIONAL_ERROR("No index specified for matrix written to 3-D or 4-D netcdf array");
    }

    int dimid[4];
    NC_CHECK(nc_inq_vardimid(ncid_, varid, dimid), varname);
    size_t len[4] = {0, 0, 0, 0};
    for (int idim = 0; idim < ndims; idim++) {
      NC_CHECK_UNNAMED(nc_inq_dimlen(ncid_, dimid[idim], len+idim));
    }
    size_t start[4] = {0, 0, 0, 0};
    if (j >= 0) { 
      start[0] = j;
      len[0] = 1;
      if (i >= 0) {
	start[1] = i;
	len[1] = 1;
      }
    }
    
    if (M.is_contiguous()) {
#ifdef REAL_IS_FLOAT
      NC_CHECK(nc_put_vara_float(ncid_, varid, start, len, M.data_pointer()),
	       varname);
#else
      NC_CHECK(nc_put_vara_double(ncid_, varid, start, len, M.data_pointer()),
	       varname);
#endif
      return true;
    }
    else {
      //      CONDITIONAL_ERROR("Cannot (yet) write non-contiguous matrices to a netcdf file");
      Matrix Mcopy;
      Mcopy.resize_contiguous(M.dimensions());
      Mcopy = M;
#ifdef REAL_IS_FLOAT
      NC_CHECK(nc_put_vara_float(ncid_, varid, start, len, Mcopy.data_pointer()),
	       varname);
#else
      NC_CHECK(nc_put_vara_double(ncid_, varid, start, len, Mcopy.data_pointer()),
	       varname);
#endif
      return true;
    }
  }
  return false;
}


bool
OutputDataFile::
write(const Matrix& M, const std::string& scope,
     const std::string& varname)
{
  if (type_ == CONFIG) {
    CONDITIONAL_ERROR("Writing matrices to cfg files is not yet implemented");
  }
  else if (type_ == NETCDF) {
    CONDITIONAL_ERROR("Cannot write a matrix to a netcdf attribute");
  }
  return false;
}


bool
OutputDataFile::
write(const Array3D& M, const std::string& varname, int j, int i)
{
  if (type_ == CONFIG) {
    if (j >= 0 || i >= 0) {
      CONDITIONAL_ERROR("Attempt to write sub-part of cfg array");
    }
    CONDITIONAL_ERROR("Writing 3D arrays to cfg files is not yet implemented");
  }
  else if (type_ == NETCDF) {
    int varid;
    NC_CHECK(nc_inq_varid(ncid_, varname.c_str(), &varid), varname);
    int ndims;
    NC_CHECK(nc_inq_varndims(ncid_, varid, &ndims), varname);
    if (ndims > 5) {
      CONDITIONAL_ERROR("Cannot write to netcdf array with more than 5 dimensions");
    }
    else if (i < 0 && ndims == 5) {
      CONDITIONAL_ERROR("No second index specified for matrix written to 5-D netcdf array");
    }
    else if (j < 0 && ndims > 3) {
      CONDITIONAL_ERROR("No index specified for matrix written to 4-D or 5-D netcdf array");
    }

    int dimid[5];
    NC_CHECK(nc_inq_vardimid(ncid_, varid, dimid), varname);
    size_t len[5] = {0, 0, 0, 0, 0};
    for (int idim = 0; idim < ndims; idim++) {
      NC_CHECK_UNNAMED(nc_inq_dimlen(ncid_, dimid[idim], len+idim));
    }
    size_t start[5] = {0, 0, 0, 0, 0};
    if (j >= 0) { 
      start[0] = j;
      len[0] = 1;
      if (i >= 0) {
	start[1] = i;
	len[1] = 1;
      }
    }
    
    if (M.is_contiguous()) {
#ifdef REAL_IS_FLOAT
      NC_CHECK(nc_put_vara_float(ncid_, varid, start, len, M.data_pointer()),
	       varname);
#else
      NC_CHECK(nc_put_vara_double(ncid_, varid, start, len, M.data_pointer()),
	       varname);
#endif
      return true;
    }
    else {
      //      CONDITIONAL_ERROR("Cannot (yet) write non-contiguous matrices to a netcdf file");
      Array3D Mcopy;
      Mcopy.resize_contiguous(M.dimensions());
      Mcopy = M;
#ifdef REAL_IS_FLOAT
      NC_CHECK(nc_put_vara_float(ncid_, varid, start, len, Mcopy.data_pointer()),
	       varname);
#else
      NC_CHECK(nc_put_vara_double(ncid_, varid, start, len, Mcopy.data_pointer()),
	       varname);
#endif
      return true;
    }
  }
  return false;
}


bool
OutputDataFile::
exist(const std::string& varname)
{
  if (type_ == CONFIG) {
    CONDITIONAL_ERROR("Unable to test for existence of variables in write-only cfg file");
  }
  else if (type_ == NETCDF) {
    int varid;
    if (nc_inq_varid(ncid_, varname.c_str(), &varid) == NC_NOERR) {
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

bool
OutputDataFile::
exist(const std::string& scope, const std::string& varname)
{
  if (type_ == CONFIG) {
    CONDITIONAL_ERROR("Unable to test for existence of variables in write-only cfg file");
  }
  else if (type_ == NETCDF) {
    int varid;
    if (nc_inq_varid(ncid_, scope.c_str(), &varid) == NC_NOERR) {
      int attid;
      if (nc_inq_attid(ncid_, varid, varname.c_str(), &attid) == NC_NOERR) {
	return true;
      }
    }
    return false;
  }
  return false;
}


// Append the current command line to the "history" attribute (NetCDF
// only)
bool
OutputDataFile::
append_history(int argc, const char* argv[])
{
  if (type_ == CONFIG) {
    CONDITIONAL_ERROR("Unable to append history in CFG file");
    return false;
  }
  else if (type_ == NETCDF) {
    size_t len;
    std::string history;
    if (nc_inq_attlen(root_ncid_, NC_GLOBAL, "history", &len) == NC_NOERR) {
      history.resize(len,' ');
      NC_CHECK(nc_get_att_text(root_ncid_, NC_GLOBAL, "history", &history[0]), 
	       "history");
      history += "\n";
    }
    std::string new_history;
    time_t t;
    std::time(&t);
    history += ctime(&t); // New-line terminated
    history[history.size()-1] = ':';
    history += " ";
    {
      // Remove directory name from command
      std::string argv0 = argv[0];
      std::string::size_type pos = argv0.find_last_of("/");
      if (pos != std::string::npos) {
	history += argv0.substr(pos+1);
      }
      else {
	history += argv0;
      }
    }

    for (int i = 1; i < argc; ++i) {
      history += " ";
      history += argv[i];
    }
    NC_CHECK(nc_put_att_text(root_ncid_, NC_GLOBAL, "history",
			     history.size(), history.c_str()), "history");
    return true;
  }
  else {
    CONDITIONAL_ERROR("Unable to append history in file of unknown type");
    return false;
  }
}


// Append the current command line to the "history" attribute (NetCDF
// only)
bool
OutputDataFile::
append_convention()
{
  if (type_ == CONFIG) {
    CONDITIONAL_ERROR("Unable to append convention in CFG file");
    return false;
  }
  else if (type_ == NETCDF) {
    size_t len;
    std::string convention;
    if (nc_inq_attlen(root_ncid_, NC_GLOBAL, "conventions", &len) == NC_NOERR) {
      convention.resize(len,' ');
      NC_CHECK(nc_get_att_text(root_ncid_, NC_GLOBAL, "conventions", &convention[0]), 
	       "conventions");
      convention += "\n";
    }
    convention += "CF-1.6";
    NC_CHECK(nc_put_att_text(root_ncid_, NC_GLOBAL, "conventions",
			     convention.size(), convention.c_str()), "conventions");
    return true;
  }
  else {
    CONDITIONAL_ERROR("Unable to append convention in file of unknown type");
    return false;
  }
}

// Append the current command line to the "history" attribute (NetCDF
// only)
bool
OutputDataFile::
append_title()
{
  if (type_ == CONFIG) {
    CONDITIONAL_ERROR("Unable to append title in CFG file");
    return false;
  }
  else if (type_ == NETCDF) {
    size_t len;
    std::string title;
    if (nc_inq_attlen(root_ncid_, NC_GLOBAL, "title", &len) == NC_NOERR) {
      title.resize(len,' ');
      NC_CHECK(nc_get_att_text(root_ncid_, NC_GLOBAL, "title", &title[0]), 
	       "title");
      title += "\n";
    }
    title += "ACM-CAP: Cloud, aerosol and precipitation from EarthCARE\'s ATLID, CPR and MSI instruments";
    NC_CHECK(nc_put_att_text(root_ncid_, NC_GLOBAL, "title",
			     title.size(), title.c_str()), "title");
    return true;
  }
  else {
    CONDITIONAL_ERROR("Unable to append title in file of unknown type");
    return false;
  }
}
