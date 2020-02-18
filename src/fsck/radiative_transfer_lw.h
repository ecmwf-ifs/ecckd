#ifndef RADIATIVE_TRANSFER_LW_H
#define RADIATIVE_TRANSFER_LW_H 1

#include <adept_arrays.h>
#include "constants.h"

/// Perform longwave radiative transfer using planck function
/// integrated across each spectral interval (dimensioned
/// half-level,wavenumber) and optical depth (dimensioned
/// layer,wavenumber) to obtain upwelling and downwelling spectral
/// fluxes (dimensioned half-level,wavenumber). The function is
/// templated so that optionally automatic differentiation can be used
/// to compute the derivative of the fluxes with respect to the
/// optical depth.
template <bool IsActive>
void
radiative_transfer_lw(const adept::Matrix& planck,                  ///< Planck function in W m-2
		      const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
		      const adept::Vector& surf_emissivity,         ///< Surface emissivity
		      const adept::Vector& surf_planck,             ///< Surface Planck function in W m-2
		      adept::Array<2,adept::Real,IsActive> flux_dn, ///< Spectral flux down in W m-2
		      adept::Array<2,adept::Real,IsActive> flux_up  ///< Spectral flux up in W m-2
		      );

/// As radiative_transfer_lw but with broadband flux outputs, to
/// reduce memory consumption if not needed
template <bool IsActive>
void
radiative_transfer_lw_bb(const adept::Matrix& planck,                  ///< Planck function in W m-2
			 const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
			 const adept::Vector& surf_emissivity,         ///< Surface emissivity
			 const adept::Vector& surf_planck,             ///< Surface Planck function in W m-2
			 adept::Array<1,adept::Real,IsActive> flux_dn, ///< Broadband flux down in W m-2
			 adept::Array<1,adept::Real,IsActive> flux_up  ///< Broadband flux up in W m-2
			 );
#endif
