#include "calc_cost_function_sw.h"
#include "radiative_transfer_sw.h"
#include "heating_rate.h"
#include "Error.h"

adept::Real
calc_cost_function_sw(adept::Real cos_sza,
		      const adept::Vector& pressure_hl,       ///< Pressure (Pa)
		      const adept::Vector& ssi,               ///< Spectral solar irradiance (W m-2)
		      adept::Real albedo,                     ///< Surface albedo
		      const adept::Matrix& bg_optical_depth,  ///< Background optical depth
		      const adept::Vector& optical_depth_fit, ///< Fitted optical depth of target gas
		      const adept::Vector& flux_dn_surf,      ///< True downwelling surface flux (W m-2)
		      const adept::Vector& flux_up_toa,       ///< True upwelling TOA flux (W m-2)
		      const adept::Matrix& hr,                ///< True heating rate (K s-1)
		      adept::Real flux_weight,                ///< Weight applied to TOA and surface fluxes
		      const adept::Vector& layer_weight,      ///< Weight applied to heating rates in each layer
		      const adept::intVector& index,          ///< Indices of wavenumbers to consider
		      bool iverbose) {
  using namespace adept;

  // Convert K s-1 to K day-1
  static const Real hr_weight = 3600.0*24.0;

  int nlay = pressure_hl.size()-1;
  int nwav = ssi.size();

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

  Vector flux_dn_fit(nlay+1);
  Vector flux_up_fit(nlay+1);
  flux_up_fit = 0.0; // By default no upwelling
  if (index.empty()) {
    if (albedo <= 0.0) {
      radiative_transfer_direct_sw_bb(cos_sza,
				      ssi,
				      bg_optical_depth,
				      optical_depth_fit,
				      flux_dn_fit);
    }
    else {
      radiative_transfer_norayleigh_sw_bb(cos_sza, ssi,
				      bg_optical_depth, optical_depth_fit, albedo,
				      flux_dn_fit, flux_up_fit);
    }
  }
  else {
    if (albedo <= 0.0) {
      radiative_transfer_direct_sw_bb(cos_sza,
				      ssi,
				      eval(bg_optical_depth(__,index)),
				      optical_depth_fit,
				      flux_dn_fit);
    }
    else {
      radiative_transfer_norayleigh_sw_bb(cos_sza, ssi,
				      eval(bg_optical_depth(__,index)), optical_depth_fit, albedo,
				      flux_dn_fit, flux_up_fit);
    }
  }

  Vector hr_fit(nlay);
 // Unallocated upwelling ensures that we don't
  // include upwelling in the heating rate calculation
  heating_rate_single(pressure_hl, flux_dn_fit, Vector(), hr_fit);

  if (iverbose) {
    adept::set_array_print_style(PRINT_STYLE_MATLAB);
    std::cerr << "    debug_partition.flux_dn_surf_true = " << flux_dn_surf_true << "\n";
    std::cerr << "    debug_partition.flux_dn_surf_fit = " << flux_dn_fit(end) << "\n";
    std::cerr << "    debug_partition.flux_up_toa_true = " << flux_up_toa_true << "\n";
    std::cerr << "    debug_partition.flux_up_toa_fit = " << flux_up_fit(0) << "\n";
    std::cerr << "    debug_partition.layer_weight = " << layer_weight << "\n";
    std::cerr << "    debug_partition.hr_true = " << hr_true << "\n";
    std::cerr << "    debug_partition.hr_fit = " << hr_fit << "\n";
    std::cerr << "    debug_partition.cf_hr = " << sqrt(hr_weight*hr_weight*sum(layer_weight*((hr_fit-hr_true)*(hr_fit-hr_true)))) << "\n";
    std::cerr << "    debug_partition.cf_flux = " << sqrt(flux_weight*((flux_dn_fit(end)-flux_dn_surf_true)*(flux_dn_fit(end)-flux_dn_surf_true)
								       +(flux_up_fit(0)-flux_up_toa_true)*(flux_up_fit(0)-flux_up_toa_true))) << "\n";
  }

  return sqrt(hr_weight*hr_weight*sum(layer_weight*((hr_fit-hr_true)*(hr_fit-hr_true)))
	      + flux_weight*((flux_dn_fit(end)-flux_dn_surf_true)*(flux_dn_fit(end)-flux_dn_surf_true)
			     +(flux_up_fit(0)-flux_up_toa_true)*(flux_up_fit(0)-flux_up_toa_true)));
}


/// Compute the cost function, in the form of the mean-squared error
/// in heating rate, of a CKD scheme
adept::aReal
calc_cost_function_ckd_sw(adept::Real cos_sza,
			  const adept::Vector& pressure_hl,       ///< Pressure (Pa)
			  const adept::Vector& ssi,               ///< Spectral solar irradiance (W m-2)
			  const adept::Vector& albedo,            ///< Spectral albedo
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
			  const adept::intVector& band_mapping,
			  adept::Vector cost_fn_per_band)
{
  using namespace adept;

  // Convert K s-1 to K day-1
  static const Real hr_weight = 3600.0*24.0;

  static bool warning_issued = false;

  int nlay = pressure_hl.size()-1;
  int ng   = optical_depth.dimension(1);

  aMatrix flux_dn_fwd_orig(nlay+1,ng);
  aMatrix flux_up_fwd_orig(nlay+1,ng);

  if (all(albedo <= 0.0)) {
    radiative_transfer_direct_sw(cos_sza, ssi,
				 optical_depth,
				 flux_dn_fwd_orig);
    flux_up_fwd_orig = 0.0;
  }
  else {
    Vector albedo_gpoint = albedo(band_mapping);
    radiative_transfer_norayleigh_sw(cos_sza, ssi,
				     optical_depth, albedo_gpoint,
				     flux_dn_fwd_orig, flux_up_fwd_orig);
  }

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
  //heating_rate(pressure_hl, flux_dn_fwd, flux_up_fwd, heating_rate_fwd);
  heating_rate(pressure_hl, flux_dn_fwd, aMatrix(), heating_rate_fwd);

  aReal cost_fn = 0.0;
  // Spectral contribution to cost function

  // Weighting for interior fluxes 
  Vector interface_weight = flux_profile_weight * 0.5*(layer_weight(range(0,end-1))+layer_weight(range(1,end)));

  Vector incoming_error_band(nband);

  //  LOG << "TOA up LBL: " << flux_up_fwd(0,__) << ", CKD: " << flux_up(0,__) << "\n";

  for (int iband = 0; iband < nband; ++iband) {
    aReal cost_fn_local;
    cost_fn_local = hr_weight*hr_weight*sum(layer_weight*(heating_rate_fwd(__,iband)-hr(__,iband))
					                *(heating_rate_fwd(__,iband)-hr(__,iband)))
      + flux_weight*((flux_dn_fwd(end,iband)-flux_dn(end,iband))*(flux_dn_fwd(end,iband)-flux_dn(end,iband))
		     +20.0*(flux_up_fwd(0,iband)-flux_up(0,iband))*(flux_up_fwd(0,iband)-flux_up(0,iband)));

    //    LOG << "  band=" << iband << ": dn cost = " << flux_weight*((flux_dn_fwd(end,iband)-flux_dn(end,iband))*(flux_dn_fwd(end,iband)-flux_dn(end,iband))) << "\n";
    //    LOG << "  band=" << iband << ": up cost = " << flux_weight*((flux_up_fwd(0,iband)-flux_up(0,iband))*(flux_up_fwd(0,iband)-flux_up(0,iband))) << "\n";

    incoming_error_band(iband) = value(flux_dn_fwd(0,iband)-flux_dn(0,iband));
    if (flux_profile_weight > 0.0) {
      cost_fn_local += sum(interface_weight*((flux_dn_fwd(range(1,end-1),iband)-flux_dn(range(1,end-1),iband))
				       *(flux_dn_fwd(range(1,end-1),iband)-flux_dn(range(1,end-1),iband))
				       +(flux_up_fwd(range(1,end-1),iband)-flux_up(range(1,end-1),iband))
				       *(flux_up_fwd(range(1,end-1),iband)-flux_up(range(1,end-1),iband))));
    }
    cost_fn += cost_fn_local;
    // Store cost function per band
    if (!cost_fn_per_band.empty()) {
      cost_fn_per_band(iband) += value(cost_fn_local);
    }
  }

  //  LOG << "BAND ERROR " << incoming_error_band << "\n";

  if (any(fabs(incoming_error_band) > 0.1) && !warning_issued) {
    WARNING << "Incoming flux in band differs by more than 0.1 W m-2: "
	    << incoming_error_band << "\n";
    ENDWARNING;
    warning_issued = true;
  }

  //  Real cost_fn_save = value(cost_fn);

  if (broadband_weight > 0.0) {
    // Broadband contribution to cost function
    cost_fn = (cost_fn*(1.0-broadband_weight))/nband
      + broadband_weight*hr_weight*hr_weight*sum(layer_weight*(sum(heating_rate_fwd-hr,1)
							      *sum(heating_rate_fwd-hr,1)));
    // Note that the sums here are over band
    cost_fn += broadband_weight*flux_weight*(sum(flux_dn_fwd(end,__)-flux_dn(end,__))
					    *sum(flux_dn_fwd(end,__)-flux_dn(end,__)));
    if (all(albedo > 0.0)) {
      cost_fn += broadband_weight*flux_weight*(sum(flux_up_fwd(0,__)-flux_up(0,__))
					      *sum(flux_up_fwd(0,__)-flux_up(0,__)));
    }
    // else:
    // Some of the bands have too much Rayleigh for the simpler RT
    // model used here to be applicable, so we can't use the broadband
    // upwelling in the cost function.
    
    if (flux_profile_weight > 0.0) {
      aVector flux_dn_error = sum(flux_dn_fwd(range(1,end-1),__)-flux_dn(range(1,end-1),__),1);
      cost_fn += broadband_weight*sum(interface_weight*(flux_dn_error*flux_dn_error));
      if (all(albedo > 0.0)) {
	aVector flux_up_error = sum(flux_up_fwd(range(1,end-1),__)-flux_up(range(1,end-1),__),1);
	cost_fn += broadband_weight*sum(interface_weight*(flux_up_error*flux_up_error));
      }
    }
  }

  if (spectral_boundary_weight > 0.0 && !spectral_flux_dn_surf.empty()) {
    cost_fn += spectral_boundary_weight
      * sum(  (flux_dn_fwd_orig(end,__)-spectral_flux_dn_surf)
	    * (flux_dn_fwd_orig(end,__)-spectral_flux_dn_surf));
  }

  return cost_fn;
}
