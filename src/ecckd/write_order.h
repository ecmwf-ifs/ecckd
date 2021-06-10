// write_order.h - Write the ordering of a spectrum to NetCDF
//
// Copyright (C) 2019- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
//
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.
//
// Author:  Robin Hogan
// Email:   r.j.hogan@ecmwf.int

#include <string>
#include <adept_arrays.h>

/// Write a NetCDF file containing the ordering of spectral intervals
/// for a single gas for use in a correlated k-distribution scheme
void
write_order(std::string& file_name,                    ///< Name of NetCDF file to write
	    int argc,                                  ///< Number of command-line args
	    const char** argv,                         ///< Command-line arguments
	    const std::string& gas_name,               ///< Formula for gas in lower case
	    const std::string& config_str,             ///< Configuration as a single string
	    const adept::Vector& band_bound1,          ///< Lower wavenumber of bands (cm-1)
	    const adept::Vector& band_bound2,          ///< Upper wavenumber of bands (cm-1)
	    const adept::Vector& wavenumber_cm_1,      ///< Wavenumber (cm-1)
	    const adept::Vector& d_wavenumber_cm_1,    ///< Wavenumber interval (cm-1)
	    const adept::intVector& iband,             ///< Band number (0 based)
	    const adept::intVector& rank,              ///< Rank of point (0 based)
	    const adept::intVector& ordered_index,     ///< Index to points in order (0 based)
	    const adept::Vector& column_optical_depth, ///< Column optical depth
	    const adept::Vector& peak_cooling_height   ///< Pseudo height of peak cooling
	    );
