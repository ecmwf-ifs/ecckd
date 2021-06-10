// read_merged_spectrum.h - Read and combine the spectral optical depths of several gases
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
#include "DataFile.h"

/// Read and combine the spectral optical depths of several gases
void
read_merged_spectrum(DataFile& config,                ///< Config file
		     int iprof,                       ///< Index of profile (0 based)
		     std::string prefix,              ///< Prefix of keys used to query "config"
		     adept::Vector& pressure_hl,      ///< Half-level pressure (Pa)
		     adept::Vector& temperature_hl,   ///< Half-level temperature (K)
		     adept::Vector& wavenumber_cm_1,  ///< Wavenumber (cm-1)
		     adept::Vector& d_wavenumber_cm_1,///< Wavenumber interval (cm-1)
		     adept::Matrix& optical_depth,    ///< Spectral optical depth 
		     std::string& molecules,          ///< Molecule formulas separated by commas
		     adept::Matrix& vmr_fl,           ///< Volume mixing ratio (mol mol-1)
		     int* ngas = 0,                   ///< Number of gases read in
		     int* ncol = 0                    ///< Number of columns in file
		     );
