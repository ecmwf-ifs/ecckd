// calc_cost_function_lw.cpp - Calculate the longwave cost function
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

#include "calc_cost_function_lw.h"
#include "radiative_transfer_lw.h"
#include "heating_rate.h"
#include "Error.h"

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
		      ) {
  using namespace adept;

  // Convert K s-1 to K day-1
  static const Real hr_weight = 3600.0*24.0;

  int nlay = pressure_hl.size()-1;
  int nwav = surf_planck.size();

  Vector hr_true;
  Real flux_dn_surf_true;
  Real flux_up_toa_true;

  if (index.empty()) {
    hr_true           = sum(hr,1);
    flux_dn_surf_true = sum(flux_dn_surf);
    flux_up_toa_true  = sum(flux_up_toa);
  }
  else {
    hr_true           = sum(hr(__,index),1);
    flux_dn_surf_true = sum(flux_dn_surf(index));
    flux_up_toa_true  = sum(flux_up_toa(index));
  }

#ifdef WASTE_MEMORY

  // I've no idea if this is still valid...

  Matrix flux_dn(nlay+1,nwav);
  Matrix flux_up(nlay+1,nwav);
  radiative_transfer_lw(planck_hl(__,index),
			eval(bg_optical_depth(__,index) + spread<1>(optical_depth_fit,nwav)),
			surf_emissivity(index),
			surf_planck(index),
			flux_dn,
			flux_up);

  Vector flux_dn_fit = sum(flux_dn,1);
  Vector flux_up_fit = sum(flux_up,1);

#else

  Vector flux_dn_fit(nlay+1);
  Vector flux_up_fit(nlay+1);
  if (index.empty()) {
    radiative_transfer_lw_bb(planck_hl,
			     bg_optical_depth,
			     optical_depth_fit,
			     surf_emissivity,
			     surf_planck,
			     flux_dn_fit,
			     flux_up_fit);
  }
  else {
    radiative_transfer_lw_bb(planck_hl(__,index),
			     eval(bg_optical_depth(__,index)),
			     optical_depth_fit,
			     surf_emissivity(index),
			     surf_planck(index),
			     flux_dn_fit,
			     flux_up_fit);
  }

#endif

  Vector hr_fit(nlay);
  heating_rate_single(pressure_hl, flux_dn_fit, flux_up_fit, hr_fit);

  //  std::cerr << "HRtrue = " << hr_true << "\n";
  //  std::cerr << "HRfit = " << hr_fit << "\n";

  return sqrt(hr_weight*hr_weight*sum(layer_weight*((hr_fit-hr_true)*(hr_fit-hr_true)))
	      + flux_weight*((flux_dn_fit(end)-flux_dn_surf_true)*(flux_dn_fit(end)-flux_dn_surf_true)
			     +(flux_up_fit(0)-flux_up_toa_true)*(flux_up_fit(0)-flux_up_toa_true)));
}


/// Compute the cost function, in the form of the mean-squared error
/// in heating rate, of a CKD scheme
adept::aReal
calc_cost_function_ckd_lw(const adept::Vector& pressure_hl,       ///< Pressure (Pa)
			  const adept::Matrix& planck_hl,         ///< Planck function (W m-2)
			  const adept::Vector& surf_emiss_orig,   ///< Surface emissivity
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
			  const adept::intVector& band_mapping)
{
  using namespace adept;

  // Convert K s-1 to K day-1
  static const Real hr_weight = 3600.0*24.0;

  int nlay = pressure_hl.size()-1;
  int ng   = optical_depth.dimension(1);

  aMatrix flux_dn_fwd_orig(nlay+1,ng);
  aMatrix flux_up_fwd_orig(nlay+1,ng);
  Vector surf_emissivity;

  if (band_mapping.empty()) {
    surf_emissivity = surf_emiss_orig;
  }
  else {
    surf_emissivity = surf_emiss_orig(band_mapping);
  }

  radiative_transfer_lw(planck_hl,
			optical_depth,
			surf_emissivity,
			surf_planck,
			flux_dn_fwd_orig,
			flux_up_fwd_orig);

  // Fluxes are to be computed relative to a reference scenario
  if (relative_ckd_flux_dn) {
    flux_dn_fwd_orig -= *relative_ckd_flux_dn;
    flux_up_fwd_orig -= *relative_ckd_flux_up;
  }


  aMatrix flux_dn_fwd, flux_up_fwd;

  int nband = ng;
  if (band_mapping.empty()) {
    flux_dn_fwd >>= flux_dn_fwd_orig;
    flux_up_fwd >>= flux_up_fwd_orig;
  }
  else {
    nband = maxval(band_mapping)+1;
    flux_dn_fwd.resize(nlay+1,nband);
    flux_up_fwd.resize(nlay+1,nband);
    for (int iband = 0; iband < nband; ++iband) {
      intVector index = find(band_mapping == iband);
      flux_dn_fwd(__,iband) = sum(flux_dn_fwd_orig(__,index),1);
      flux_up_fwd(__,iband) = sum(flux_up_fwd_orig(__,index),1);
    }
  }

  aMatrix heating_rate_fwd(nlay,nband);
  heating_rate(pressure_hl, flux_dn_fwd, flux_up_fwd, heating_rate_fwd);

  aReal cost_fn = 0.0;
  // Spectral contribution to cost function

  // Weighting for interior fluxes 
  Vector interface_weight = flux_profile_weight * 0.5*(layer_weight(range(0,end-1))+layer_weight(range(1,end)));

  for (int iband = 0; iband < nband; ++iband) {
    cost_fn += hr_weight*hr_weight*sum(layer_weight*((heating_rate_fwd(__,iband))-hr(__,iband))
				       *(heating_rate_fwd(__,iband)-hr(__,iband)))
      + flux_weight*((flux_dn_fwd(end,iband)-flux_dn(end,iband))*(flux_dn_fwd(end,iband)-flux_dn(end,iband))
		     +(flux_up_fwd(0,iband)-flux_up(0,iband))*(flux_up_fwd(0,iband)-flux_up(0,iband)));
    if (flux_profile_weight > 0.0) {
      cost_fn += sum(interface_weight*((flux_dn_fwd(range(1,end-1),iband)-flux_dn(range(1,end-1),iband))
				       *(flux_dn_fwd(range(1,end-1),iband)-flux_dn(range(1,end-1),iband))
				       +(flux_up_fwd(range(1,end-1),iband)-flux_up(range(1,end-1),iband))
				       *(flux_up_fwd(range(1,end-1),iband)-flux_up(range(1,end-1),iband))));
    }
  }

  // Broadband contribution to cost function
  cost_fn = (cost_fn*(1.0-broadband_weight))/nband
    + broadband_weight*hr_weight*hr_weight*sum(layer_weight*(sum(heating_rate_fwd-hr,1)
							     *sum(heating_rate_fwd-hr,1)))
    + broadband_weight*flux_weight*(sum(flux_dn_fwd(end,__)-flux_dn(end,__))*sum(flux_dn_fwd(end,__)-flux_dn(end,__))
				    +sum(flux_up_fwd(0,__)-flux_up(0,__))*sum(flux_up_fwd(0,__)-flux_up(0,__)));
  //LOG << sum(flux_dn_fwd(end,__)) << " " << sum(flux_dn_fwd_orig(end,__)) << " " << sum(flux_dn(end,__)) << "\n";

  if (flux_profile_weight > 0.0) {
    aVector flux_dn_error = sum(flux_dn_fwd(range(1,end-1),__)-flux_dn(range(1,end-1),__),1);
    aVector flux_up_error = sum(flux_up_fwd(range(1,end-1),__)-flux_up(range(1,end-1),__),1);
    cost_fn += broadband_weight*sum(interface_weight*(flux_dn_error*flux_dn_error
						      +flux_up_error*flux_up_error));
  }

  if (spectral_boundary_weight > 0.0 && !spectral_flux_dn_surf.empty() && !spectral_flux_up_toa.empty()) {
    cost_fn += spectral_boundary_weight
      * sum(  (flux_dn_fwd_orig(end,__)-spectral_flux_dn_surf)
	     *(flux_dn_fwd_orig(end,__)-spectral_flux_dn_surf)
	    + (flux_up_fwd_orig(0  ,__)-spectral_flux_up_toa )
	     *(flux_up_fwd_orig(0  ,__)-spectral_flux_up_toa ));
  }

  return cost_fn;
}
