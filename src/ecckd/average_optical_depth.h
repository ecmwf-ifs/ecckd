// average_optical_depth.h - Average optical depth in a spectral interval
//
// Copyright (C) 2020- ECMWF.
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

#ifndef AVERAGE_OPTICAL_DEPTH_H
#define AVERAGE_OPTICAL_DEPTH_H

// Average optical depths for each pressure level to each g point
void
average_optical_depth_to_g_point(int ng,                             ///< Number of g points
				 adept::Real reference_surface_vmr,  ///< Volume mixing ratio (mol mol-1)
				 const adept::Vector& pressure_fl,   ///< Full-level pressure (Pa)
				 const adept::Vector& pressure_hl,   ///< Half-level pressure (Pa)
				 const adept::Vector& g_point,       ///< G point for each wavenumber
				 const adept::Matrix& optical_depth, ///< Optical depth (pressure,wavenumber)
				 const adept::Matrix& planck_fl,     ///< Planck function, W m-2 (pressure,wavenumber)
				 const std::string& averaging_method,
				 adept::Matrix molar_abs);           ///< Molar absorption coefficient, m2 mol-1 (pressure,g-point)
#endif
