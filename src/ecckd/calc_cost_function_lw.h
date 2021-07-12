// calc_cost_function_lw.h - Calculate the longwave cost function
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

#ifndef CALC_COST_FUNCTION_LW_H
#define CALC_COST_FUNCTION_LW_H 1

#include <adept_arrays.h>

/// Compute the cost function, in the form of the root-mean-squared
/// error in heating rate, associated with averaging the optical depth
/// of the target gas across a range of wavenumbers
adept::Real
calc_cost_function_lw(const adept::Vector& pressure_hl,       ///< Pressure (Pa)
		      const adept::Matrix& planck_hl,         ///< Planck function (W m-2)
		      const adept::Vector& surf_emissivity,   ///< Surface emissivity
		      const adept::Vector& surf_planck,       ///< Surface spectral Planck function (W m-2)
		      const adept::Matrix& bg_optical_depth,  ///< Background optical depth
		      const adept::Vector& optical_depth_fit, ///< Fitted optical depth of target gas
		      const adept::Vector& flux_dn_surf,      ///< True downwelling surface flux (W m-2)
		      const adept::Vector& flux_up_toa,       ///< True upwelling TOA flux (W m-2)
		      const adept::Matrix& hr,                ///< True heating rate (K s-1)
		      adept::Real flux_weight,                ///< Weight applied to TOA and surface fluxes
		      const adept::Vector& layer_weight,      ///< Weight applied to heating rates in each layer
		      const adept::intVector& index = adept::intVector() ///< Indices of wavenumbers to consider
		      );

/// Compute the cost function, in the form of the mean-squared error
/// in heating rate, of a CKD scheme
adept::aReal
calc_cost_function_ckd_lw(const adept::Vector& pressure_hl,       ///< Pressure (Pa)
			  const adept::Matrix& planck_hl,         ///< Planck function (W m-2)
			  const adept::Vector& surf_emissivity,   ///< Surface emissivity
			  const adept::Vector& surf_planck,       ///< Surface spectral Planck function (W m-2)
			  const adept::aMatrix& optical_depth,    ///< Optical depth of gases
			  const adept::Matrix& flux_dn,           ///< True downwelling flux (W m-2)
			  const adept::Matrix& flux_up,           ///< True upwelling flux (W m-2)
			  const adept::Matrix& hr,                ///< True heating rate (K s-1)
			  const adept::Vector& spectral_flux_dn_surf, ///< g-point surface downward flux (W m-2)
			  const adept::Vector& spectral_flux_up_toa,  ///< g-point TOA upward flux (W m-2)
			  adept::Real flux_weight,                ///< Weight applied to TOA and surface fluxes
			  adept::Real flux_profile_weight,        ///< Weight applied to other fluxes
			  adept::Real broadband_weight,           ///< Weight of broadband vs spectral (0-1)
			  adept::Real spectral_boundary_weight,   ///< Weight of spectral boundary fluxes
			  const adept::Vector& layer_weight,      ///< Weight applied to heating rates in each layer
			  adept::Matrix* relative_ckd_flux_dn,    ///< Subtract relative-to flux dn, if not NULL
			  adept::Matrix* relative_ckd_flux_up,    ///< Subtract relative-to flux up, if not NULL
			  const adept::intVector& band_mapping = adept::intVector());

#endif
