#ifndef LBL_FLUXES_H
#define LBL_FLUXES_H 1

#include <string>
#include <vector>
#include <adept_arrays.h>

struct LblFluxes {

  LblFluxes() { }

  LblFluxes(const std::string& file_name,
	    const adept::intVector& band_mapping = adept::intVector(), ///< Mapping from g-point to band
	    const adept::intVector& g_point = adept::intVector()) { ///< Mapping from wavenumber to g-point
    read(file_name, band_mapping, g_point);
  }

  adept::Matrix pressure_hl_; ///< Pressure (Pa), dimensioned (column, half_level)
  adept::Matrix temperature_hl_; ///< Temperature (K), dimensioned (column, half_level)
  adept::Array3D vmr_fl_; ///< Volume mixing ratio (mol mol-1), dimensioned (column, gas, level)
  adept::Matrix flux_up_; ///< Upwelling flux (W m-2), dimensioned (column, half_level)
  adept::Matrix flux_dn_; ///< Downwelling flux (W m-2), dimensioned (column, half_level)
  adept::Array3D spectral_flux_up_; ///< Upwelling flux (W m-2), dimensioned (column, half_level, g_point/band)
  adept::Array3D spectral_flux_dn_; ///< Downwelling flux (W m-2), dimensioned (column, half_level, g_point/band)
  adept::Matrix spectral_flux_dn_surf_; ///< Downwelling surface flux (W m-2) at g-points, dimensioned (column, g_point)
  adept::Matrix spectral_flux_up_toa_; ///< Upwelling TOA flux (W m-2) at g-points, dimensioned (column, g_point)
  //  adept::Array3D band_flux_up_; ///< Upwelling flux (W m-2), dimensioned (column, half_level, band)
  //  adept::Array3D band_flux_dn_; ///< Downwelling flux (W m-2), dimensioned (column, half_level, band)
  adept::Matrix heating_rate_; /// Heating rate (K d-1), dimensioned (column, level)
  adept::Array3D spectral_heating_rate_; /// Spectral heating rate (K d-1), dimensioned (column, level, g_point)
  //  adept::Array3D band_heating_rate_; /// Spectral heating rate (K d-1), dimensioned (column, level, band)
  adept::Vector mu0_; ///< Cosine of solar zenith angle
  adept::Vector effective_spectral_albedo_; // Surface albedo in each band, actually up/direct_dn

  adept::Matrix surf_emissivity_, surf_planck_;
  adept::Array3D planck_hl_;
  adept::Vector solar_irradiance_; // At g points

  std::vector<std::string> molecules_;

  // Total solar irradiance (W m-2)
  adept::Real tsi_;

  int ngas() { return molecules_.size(); }

  int nspec() { return spectral_flux_up_.size(2); }

  void read(const std::string& file_name,
	    const adept::intVector& band_mapping = adept::intVector(),
	    const adept::intVector& g_point = adept::intVector());

  // For optimizing the coefficients of minor greenhouse gases, we
  // need to compute forcing relative to a reference, so need the
  // ability to subtract one set of fluxes and heating rates from
  // another
  void subtract(const LblFluxes& source);

  // Fill the gas_mapping vector to map from the CKD gas indices to
  // the LBL concentrations
  void make_gas_mapping(const std::vector<std::string>& molecules);

  /// Calculate fluxes given a 3D array of optical depths for a CKD
  /// model for the same scenario as the present object
  void calc_ckd_fluxes(const adept::Array3D& optical_depth,
		       adept::Array3D& flux_dn, adept::Array3D& flux_up) const;

  // Remove upwelling fluxes affected by Rayleigh scattering since
  // they will not be adequately modelled by the fast scheme used for
  // the optimization
  void mask_rayleigh_up(adept::Real max_no_rayleigh_wavenumber);

  // Do we have fluxes in either g-points or bands?
  bool have_spectral_fluxes = false;

  // Are those fluxes in bands (rather than g-points)?
  bool have_band_fluxes = false;

  // Lower and upper wavenumbers (cm-1) for bands 
  adept::Vector band_wavenumber1_, band_wavenumber2_;

  // Band index for each g point
  adept::intVector iband_per_g;

  // Mapping between gas index of CKD model and gas index of
  // mole_fraction_fl in LBL file (stored as vmr_fl_ in this struct),
  // with -1 indicating a missing gas (e.g. "composite" which
  // corresponds to a gas mixture)
  adept::intVector gas_mapping;

  bool is_sw() const { return is_sw_; }

  bool is_sw_;

};


#endif
