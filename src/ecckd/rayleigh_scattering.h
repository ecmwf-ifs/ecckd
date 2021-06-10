// rayleigh_scattering.h - Compute Rayleigh molar scattering coefficient
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

#ifndef RAYLEIGH_SCATTERING_H
#define RAYLEIGH_SCATTERING_H 1

#include <adept_arrays.h>
#include "constants.h"

/// Compute the Rayleigh molar scattering coefficient (m2 mol-1) as a
/// function of wavenumber (cm-1)
inline
adept::Vector
rayleigh_molar_scattering_coeff(const adept::Vector& wavenumber_cm_1) {

  using namespace adept;

  // Convert input to wavelength in microns
  Vector wavelength_um = 10000.0 / wavenumber_cm_1;
  Vector scat_coeff(wavenumber_cm_1.size());
  // Use Bucholtz (1995) model to get per-molecule cross-section in m2
  scat_coeff.where(wavelength_um < 0.5)
    = either_or(3.01577e-32 * pow(wavelength_um, -( 3.55212
						    + 1.35579 * wavelength_um
						    + 0.11563 / wavelength_um)),
		4.01061e-32 * pow(wavelength_um, -( 3.99668
						    + 0.00110298 * wavelength_um
						    + 0.0271393  / wavelength_um)));
  // Convert to molar cross-section
  scat_coeff *= AVOGADRO_CONSTANT;
  return scat_coeff;
}


#endif
