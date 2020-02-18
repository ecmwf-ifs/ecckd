#ifndef CALC_COST_FUNCTION_LW_H
#define CALC_COST_FUNCTION_LW_H 1

#include <adept_arrays.h>
#include "radiative_transfer_lw.h"
#include "heating_rate.h"


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
		      const adept::intVector& index           ///< Indices of wavenumbers to consider
		      );

/// Compute the cost function, in the form of the mean-squared error
/// in heating rate, of a CKD scheme
adept::aReal
calc_cost_function_lw(const adept::Vector& pressure_hl,       ///< Pressure (Pa)
		      const adept::Matrix& planck_hl,         ///< Planck function (W m-2)
		      const adept::Vector& surf_emissivity,   ///< Surface emissivity
		      const adept::Vector& surf_planck,       ///< Surface spectral Planck function (W m-2)
		      const adept::aMatrix& optical_depth,    ///< Optical depth of gases
		      const adept::Matrix& flux_dn,           ///< True downwelling flux (W m-2)
		      const adept::Matrix& flux_up,           ///< True upwelling flux (W m-2)
		      const adept::Matrix& hr,                ///< True heating rate (K s-1)
		      adept::Real flux_weight,                ///< Weight applied to TOA and surface fluxes
		      adept::Real flux_profile_weight,        ///< Weight applied to other fluxes
		      adept::Real broadband_weight,           ///< Weight of broadband vs spectral (0-1)
		      const adept::Vector& layer_weight       ///< Weight applied to heating rates in each layer
		      );

#endif
