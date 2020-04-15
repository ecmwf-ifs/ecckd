#include "lbl_fluxes.h"
#include "heating_rate.h"
#include "planck_function.h"
#include "DataFile.h"

using namespace adept;

void
LblFluxes::read(const std::string& file_name, const intVector& band_mapping)
{
  LOG << "Reading " << file_name << "\n";
  DataFile file(file_name);
  file.read(pressure_hl_, "pressure_hl");
  file.read(temperature_hl_, "temperature_hl");
  file.read(vmr_fl_, "mole_fraction_fl");

  file.read(flux_dn_, "flux_dn_lw");
  file.read(flux_up_, "flux_up_lw");

  if (file.exist("spectral_flux_up_lw")) {
    file.read(spectral_flux_up_, "spectral_flux_up_lw");
    file.read(spectral_flux_dn_, "spectral_flux_dn_lw");
    have_spectral_fluxes = true;
    have_band_fluxes = false;
  }
  else if (file.exist("band_flux_up_lw")) {
    if (band_mapping.empty()) {
      file.read(spectral_flux_up_, "band_flux_up_lw");
      file.read(spectral_flux_dn_, "band_flux_dn_lw");
      file.read(band_wavenumber1_, "band_wavenumber1_lw");
      file.read(band_wavenumber2_, "band_wavenumber2_lw");
    }
    else {
      // Optionally map from narrow to wide bands
      Array3D spectral_flux;
      int nband = maxval(band_mapping)+1;
      file.read(spectral_flux, "band_flux_up_lw");
      LOG << "  Mapping fluxes from " << spectral_flux.size(2) 
	  << " to " << nband << " bands\n";
      spectral_flux_up_.resize(spectral_flux.size(0), spectral_flux.size(1), nband);
      for (int jband = 0; jband < nband; ++jband) {
	spectral_flux_up_(__,__,jband) = sum(spectral_flux(__,__,find(band_mapping==jband)),2);
      }
      file.read(spectral_flux, "band_flux_dn_lw");
      spectral_flux_dn_.resize(spectral_flux.size(0), spectral_flux.size(1), nband);
      for (int jband = 0; jband < nband; ++jband) {
	spectral_flux_dn_(__,__,jband) = sum(spectral_flux(__,__,find(band_mapping==jband)),2);
      }
      Vector band_wn1, band_wn2;
      band_wavenumber1_.resize(nband);
      band_wavenumber2_.resize(nband);
      file.read(band_wn1, "band_wavenumber1_lw");
      file.read(band_wn2, "band_wavenumber2_lw");
      for (int jband = 0; jband < nband; ++jband) {
	band_wavenumber1_(jband) = minval(band_wn1(find(band_mapping==jband)));
	band_wavenumber2_(jband) = maxval(band_wn2(find(band_mapping==jband)));
      }
    }
    have_spectral_fluxes = true;
    have_band_fluxes = true;
  }

  std::string molecules_str, molecule;
  file.read(molecules_str, DATA_FILE_GLOBAL_SCOPE, "constituent_id");
  std::stringstream molecules_s(molecules_str);
  while (std::getline(molecules_s, molecule, ' ')) {
    molecules_.push_back(molecule);
  }

  int ncol = pressure_hl_.dimension(0);
  int nlev = pressure_hl_.dimension(1)-1;
  int nspec   = spectral_flux_up_.dimension(2);
  //  int nband= band_flux_up_.dimension(2);

  heating_rate_.resize(ncol,nlev);
  for (int icol = 0; icol < ncol; ++icol) {
    heating_rate_single(pressure_hl_(icol,__), flux_dn_(icol,__), flux_up_(icol,__),
		 heating_rate_(icol,__));
  }

  if (have_spectral_fluxes) {
    spectral_heating_rate_.resize(ncol,nlev,nspec);
    for (int icol = 0; icol < ncol; ++icol) {
      heating_rate(pressure_hl_(icol,__), spectral_flux_dn_(icol,__,__), spectral_flux_up_(icol,__,__),
		   spectral_heating_rate_(icol,__,__));
    }
  }

  /*
  if (have_band_fluxes) {
    band_heating_rate_.resize(ncol,nlev,nspec);
    for (int icol = 0; icol < ncol; ++icol) {
      heating_rate(pressure_hl_(icol,__), band_flux_dn_(icol,__,__), band_flux_up_(icol,__,__),
		   band_heating_rate_(icol,__,__));
    }
  }
  */

  surf_emissivity_.resize(ncol,nspec);
  surf_emissivity_ = 1.0;

}
