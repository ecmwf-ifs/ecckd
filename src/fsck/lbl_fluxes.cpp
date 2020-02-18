#include "lbl_fluxes.h"
#include "heating_rate.h"
#include "planck_function.h"
#include "DataFile.h"

using namespace adept;

void
LblFluxes::read(const std::string& file_name)
{
  DataFile file(file_name);
  file.read(pressure_hl_, "pressure_hl");
  file.read(temperature_hl_, "temperature_hl");
  file.read(vmr_fl_, "vmr_fl");
  file.read(spectral_flux_up_, "spectral_flux_up_lw");
  file.read(spectral_flux_dn_, "spectral_flux_dn_lw");

  std::string molecules_str, molecule;
  file.read(molecules_str, DATA_FILE_GLOBAL_SCOPE, "molecules");
  std::stringstream molecules_s(molecules_str);
  while (std::getline(molecules_s, molecule, ' ')) {
    molecules_.push_back(molecule);
  }

  flux_dn_ = sum(spectral_flux_dn_,2);
  flux_up_ = sum(spectral_flux_up_,2);

  int ncol = pressure_hl_.dimension(0);
  int nlev = pressure_hl_.dimension(1)-1;
  int ng   = spectral_flux_up_.dimension(2);
  spectral_heating_rate_.resize(ncol,nlev,ng);

  for (int icol = 0; icol < ncol; ++icol) {
    heating_rate(pressure_hl_(icol,__), spectral_flux_dn_(icol,__,__), spectral_flux_up_(icol,__,__),
		 spectral_heating_rate_(icol,__,__));
  }
  heating_rate_ = sum(spectral_heating_rate_,2);

  surf_emissivity_.resize(ncol,ng);
  surf_emissivity_ = 1.0;

}
