// planck_function.cpp - Compute Planck function
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

#include "planck_function.h"
#include "Error.h"

/// Compute the Planck function for each temperature and wavenumber as
/// the integral across a wavenumber interval, as a spectral
/// irradiance in W m-2
void
planck_function(const adept::Vector& temperature,       ///< Temperature in K
		const adept::Vector& wavenumber_cm_1,   ///< Wavenumber in cm-1
		const adept::Vector& d_wavenumber_cm_1, ///< Wavenumber interval in cm-1
		adept::Matrix ans                       ///< Planck function in W m-2
		) {
  using namespace adept;

  static const Real h = 6.62606896e-34;
  static const Real c = 2.99792458e8;
  static const Real k = 1.3806504e-23;
  static const Real inv_cm_2_Hz = 100.0*c;
  static const Real pi = 3.14159265358979323846;

  if (wavenumber_cm_1.size() != d_wavenumber_cm_1.size()) {
    ERROR << "wavenumber_cm_1 and d_wavenumber_cm_1 must have the same length";
    THROW(PARAMETER_ERROR);
  }

  if (ans.size(0) != temperature.size() || ans.size(1) != wavenumber_cm_1.size()) {
    ERROR << "planck_function output array is " << ans.size(0) << "x" << ans.size(1)
	  << " but should be " << temperature.size() << "x" << wavenumber_cm_1.size();
    THROW(PARAMETER_ERROR);
  }

  //  ans.resize(temperature.size(), wavenumber_cm_1.size());
  Vector freq = wavenumber_cm_1*inv_cm_2_Hz;
  Vector prefactor = (d_wavenumber_cm_1*2.0*h*inv_cm_2_Hz*pi/(c*c))
		      * (freq*freq*freq);
#pragma omp parallel for
  for (int i = 0; i < temperature.size(); i++) {
    ans.soft_link()(i,__) = prefactor / (exp((h/k)*(freq/temperature(i)))-1.0);
  }
}
