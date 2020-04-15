#ifndef LBL_FLUXES_H
#define LBL_FLUXES_H 1

#include <string>
#include <vector>
#include <adept_arrays.h>

struct LblFluxes {

  LblFluxes(const std::string& file_name, const adept::intVector& band_mapping = adept::intVector()) {
    read(file_name, band_mapping);
  }

  adept::Matrix pressure_hl_; ///< Pressure (Pa), dimensioned (column, half_level)
  adept::Matrix temperature_hl_; ///< Temperature (K), dimensioned (column, half_level)
  adept::Array3D vmr_fl_; ///< Volume mixing ratio (mol mol-1), dimensioned (column, gas, level)
  adept::Matrix flux_up_; ///< Upwelling flux (W m-2), dimensioned (column, half_level)
  adept::Matrix flux_dn_; ///< Downwelling flux (W m-2), dimensioned (column, half_level)
  adept::Array3D spectral_flux_up_; ///< Upwelling flux (W m-2), dimensioned (column, half_level, g_point/band)
  adept::Array3D spectral_flux_dn_; ///< Downwelling flux (W m-2), dimensioned (column, half_level, g_point/band)
  //  adept::Array3D band_flux_up_; ///< Upwelling flux (W m-2), dimensioned (column, half_level, band)
  //  adept::Array3D band_flux_dn_; ///< Downwelling flux (W m-2), dimensioned (column, half_level, band)
  adept::Matrix heating_rate_; /// Heating rate (K d-1), dimensioned (column, level)
  adept::Array3D spectral_heating_rate_; /// Spectral heating rate (K d-1), dimensioned (column, level, g_point)
  //  adept::Array3D band_heating_rate_; /// Spectral heating rate (K d-1), dimensioned (column, level, band)

  adept::Matrix surf_emissivity_, surf_planck_;
  adept::Array3D planck_hl_;

  std::vector<std::string> molecules_;

  int ngas() { return molecules_.size(); }

  int nspec() { return spectral_flux_up_.size(2); }

  void read(const std::string& file_name, const adept::intVector& band_mapping = adept::intVector());

  // Do we have fluxes in either g-points or bands?
  bool have_spectral_fluxes = false;

  // Are those fluxes in bands (rather than g-points)?
  bool have_band_fluxes = false;

  // Lower and upper wavenumbers (cm-1) for bands 
  adept::Vector band_wavenumber1_, band_wavenumber2_;

  // Band index for each g point
  adept::intVector iband_per_g;

};


#endif
