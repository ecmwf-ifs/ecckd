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
  int nwav = index.size();

  Vector hr_true         = sum(hr(__,index),1);
  Real flux_dn_surf_true = sum(flux_dn_surf(index));
  Real flux_up_toa_true  = sum(flux_up_toa(index));

#ifdef WASTE_MEMORY

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
  radiative_transfer_lw_bb(planck_hl(__,index),
			   eval(bg_optical_depth(__,index) + spread<1>(optical_depth_fit,nwav)),
			   surf_emissivity(index),
			   surf_planck(index),
			   flux_dn_fit,
			   flux_up_fit);

#endif

  Vector hr_fit(nlay);
  heating_rate_single(pressure_hl, flux_dn_fit, flux_up_fit, hr_fit);

  return sqrt(hr_weight*hr_weight*sum(layer_weight*((hr_fit-hr_true)*(hr_fit-hr_true)))
	      + flux_weight*((flux_dn_fit(end)-flux_dn_surf_true)*(flux_dn_fit(end)-flux_dn_surf_true)
			     +(flux_up_fit(0)-flux_up_toa_true)*(flux_up_fit(0)-flux_up_toa_true)));
}


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
		      )
{
  using namespace adept;

  // Convert K s-1 to K day-1
  static const Real hr_weight = 3600.0*24.0;

  int nlay = pressure_hl.size()-1;
  int ng   = optical_depth.dimension(1);

  aMatrix flux_dn_fwd(nlay+1,ng);
  aMatrix flux_up_fwd(nlay+1,ng);
  radiative_transfer_lw(planck_hl,
			optical_depth,
			surf_emissivity,
			surf_planck,
			flux_dn_fwd,
			flux_up_fwd);
  aMatrix heating_rate_fwd(nlay,ng);
  heating_rate(pressure_hl, flux_dn_fwd, flux_up_fwd, heating_rate_fwd);

  aReal cost_fn = 0.0;
  // Spectral contribution to cost function

  /*
  LOG << "hr = " << hr << "\n";
  LOG << "flux_dn_surf = " << flux_dn_surf << "\n";
  LOG << "flux_up_toa  = " << flux_up_toa << "\n";
  LOG << "heating_rate_fwd = " << heating_rate_fwd << "\n";
  LOG << "flux_dn_fwd = " << flux_dn_fwd << "\n";
  LOG << "flux_up_fwd = " << flux_up_fwd << "\n";
  */

  // Weighting for interior fluxes 
  Vector interface_weight = flux_profile_weight * 0.5*(layer_weight(range(0,end-1))+layer_weight(range(1,end)));

  for (int ig = 0; ig < ng; ++ig) {
    cost_fn += hr_weight*hr_weight*sum(layer_weight*((heating_rate_fwd(__,ig))-hr(__,ig))
				                    *(heating_rate_fwd(__,ig)-hr(__,ig)))
      + flux_weight*((flux_dn_fwd(end,ig)-flux_dn(end,ig))*(flux_dn_fwd(end,ig)-flux_dn(end,ig))
		     +(flux_up_fwd(0,ig)-flux_up(0,ig))*(flux_up_fwd(0,ig)-flux_up(0,ig)));
    if (flux_profile_weight > 0.0) {
      cost_fn += sum(interface_weight*((flux_dn_fwd(range(1,end-1),ig)-flux_dn(range(1,end-1),ig))
				       *(flux_dn_fwd(range(1,end-1),ig)-flux_dn(range(1,end-1),ig))
				       +(flux_up_fwd(range(1,end-1),ig)-flux_up(range(1,end-1),ig))
				       *(flux_up_fwd(range(1,end-1),ig)-flux_up(range(1,end-1),ig))));
    }
  }
  // Broadband contribution to cost function
  cost_fn = (cost_fn*(1.0-broadband_weight))/ng
    + broadband_weight*hr_weight*hr_weight*sum(layer_weight*(sum(heating_rate_fwd-hr,1)
							     *sum(heating_rate_fwd-hr,1)))
    + broadband_weight*flux_weight*(sum(flux_dn_fwd(end,__)-flux_dn(end,__))*sum(flux_dn_fwd(end,__)-flux_dn(end,__))
				    +sum(flux_up_fwd(0,__)-flux_up(0,__))*sum(flux_up_fwd(0,__)-flux_up(0,__)));
  if (flux_profile_weight > 0.0) {
    aVector flux_dn_error = sum(flux_dn_fwd(range(1,end-1),__)-flux_dn(range(1,end-1),__),1);
    aVector flux_up_error = sum(flux_up_fwd(range(1,end-1),__)-flux_up(range(1,end-1),__),1);
    cost_fn += broadband_weight*sum(interface_weight*(flux_dn_error*flux_dn_error
						      +flux_up_error*flux_up_error));
  }
  return cost_fn;
}
