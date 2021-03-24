#ifndef RADIATIVE_TRANSFER_SW_H
#define RADIATIVE_TRANSFER_SW_H 1

#include <adept_arrays.h>
#include "constants.h"

/// Perform direct-beam shortwave radiative transfer using planck function
/// integrated across each spectral interval (dimensioned
/// half-level,wavenumber) and optical depth (dimensioned
/// layer,wavenumber) to obtain upwelling and downwelling spectral
/// fluxes (dimensioned half-level,wavenumber). The function is
/// templated so that optionally automatic differentiation can be used
/// to compute the derivative of the fluxes with respect to the
/// optical depth.
template <bool IsActive>
void
radiative_transfer_direct_sw(adept::Real cos_sza,             ///< Cosine of the solar zenith angle
			     const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
			     const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
			     adept::Array<2,adept::Real,IsActive> flux_dn ///< Spectral flux down in W m-2
			     );

/// As radiative_transfer_sw but with broadband flux outputs, to
/// reduce memory consumption if not needed
template <bool IsActive>
void
radiative_transfer_direct_sw_bb(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
			 const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
			 const adept::Array<2,adept::Real,IsActive>& spectral_od, ///< Spectral layer optical depth
			 const adept::Array<1,adept::Real,IsActive>& grey_od,     ///< Grey layer optical depth
			 adept::Array<1,adept::Real,IsActive> flux_dn ///< Broadband flux down in W m-2
			 );

/// As radiative_transfer_direct_sw but the upwelling flux is also
/// computed assuming no Rayleigh scattering
template <bool IsActive>
void
radiative_transfer_norayleigh_sw(adept::Real cos_sza,      ///< Cosine of the solar zenith angle
				 const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
				 const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
				 const adept::Vector& albedo, ///< Surface albedo
				 adept::Array<2,adept::Real,IsActive> flux_dn, ///< Spectral flux down in W m-2
				 adept::Array<2,adept::Real,IsActive> flux_up ///< Spectral flux up in W m-2
				 );

/// As radiative_transfer_norayleigh_sw but with broadband outputs
template <bool IsActive>
void
radiative_transfer_norayleigh_sw_bb(adept::Real cos_sza,      ///< Cosine of the solar zenith angle
				 const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
				 const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
				 const adept::Array<1,adept::Real,IsActive>& grey_od,     ///< Grey layer optical depth
				 adept::Real albedo, ///< Surface albedo
				 adept::Array<1,adept::Real,IsActive> flux_dn, ///< Spectral flux down in W m-2
				 adept::Array<1,adept::Real,IsActive> flux_up ///< Spectral flux up in W m-2
				 );

#endif
