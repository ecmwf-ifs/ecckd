// lbl_fluxes.cpp - Manipulate line-by-line fluxes
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

#include "lbl_fluxes.h"
#include "heating_rate.h"
#include "planck_function.h"
#include "DataFile.h"
#include "radiative_transfer_lw.h"
#include "radiative_transfer_sw.h"

using namespace adept;

static
void
repeat_matrix(Matrix& mat, int nrep)
{
  int ncol = mat.size(0);
  int nlev = mat.size(1);
  Matrix newmat(ncol*nrep,nlev);
  for (int icol = 0; icol < ncol; ++icol) {
    newmat(range(icol*nrep,(icol+1)*nrep-1),__) = spread<0>(mat(icol,__),nrep);
  }
  swap(mat,newmat);
}

static
void
repeat_array3D(Array3D& mat, int nrep)
{
  int ncol = mat.size(0);
  int ngas = mat.size(1);
  int nlev = mat.size(2);
  Array3D newmat(ncol*nrep,ngas,nlev);
  for (int icol = 0; icol < ncol; ++icol) {
    newmat(range(icol*nrep,(icol+1)*nrep-1),__,__) = spread<0>(mat(icol,__,__),nrep);
  }
  swap(mat,newmat);
}

void
LblFluxes::read(const std::string& file_name, const intVector& band_mapping,
		const adept::intVector& g_point)
{
  LOG << "Reading LBL fluxes from " << file_name << "\n";

  is_sw_ = false;
  std::string domain_str = "lw";

  DataFile file(file_name);
  file.read(pressure_hl_, "pressure_hl");
  file.read(temperature_hl_, "temperature_hl");
  file.read(vmr_fl_, "mole_fraction_fl");

  int ncol = pressure_hl_.dimension(0);
  int nlev = pressure_hl_.dimension(1)-1;

  if (file.exist("mu0")) {
    is_sw_ = true;
    domain_str = "sw";
  }

  if (is_sw_) {
    // Read cosine of solar zenith angle
    Vector mu0_all;
    file.read(mu0_all, "mu0");

    // We can have all sun angles...
    //    intVector index_sza = range(0,mu0_all.size()-1);
    // ...or just 60 degrees
    //    intVector index_sza = {2};
    intVector index_sza = {0, 2, 4};

    // Number of sun angles to consider
    int nsza = index_sza.size();

    int ncol_new = ncol*nsza;

    if (nsza > 1) {
      // Need to spread the dimensions of the pressure, temperature
      // and volume mixing ratio
      repeat_matrix(pressure_hl_, nsza);
      repeat_matrix(temperature_hl_, nsza);
      repeat_array3D(vmr_fl_, nsza);
    }

    Array3D flux_up_tmp, flux_dn_tmp;
    file.read(flux_dn_tmp, "flux_dn_direct_sw");
    file.read(flux_up_tmp, "flux_up_sw");

    // Pack cosine of solar zenith angle into a vector
    mu0_.resize(ncol_new);
    mu0_.reshape(ncol,nsza) = spread<0>(mu0_all(index_sza),ncol);

    // Pack fluxes from each column and angle into one dimension
    flux_dn_.resize(ncol_new,nlev+1);
    flux_up_.resize(flux_dn_.dimensions());

    int icol_new = 0;
    for (int icol = 0; icol < ncol; ++icol) {
      for (int isza = 0; isza < nsza; ++isza) {
	flux_dn_(icol_new,__) = flux_dn_tmp(icol,index_sza(isza),__);
	flux_up_(icol_new++,__) = flux_up_tmp(icol,index_sza(isza),__);
      }
    }

    // Save total solar irradiance
    tsi_ = flux_dn_(0,0) / mu0_(0);

    Array4D spectral_flux_up, spectral_flux_dn;
    if (file.exist("spectral_flux_dn_direct_sw")) {
      file.read(spectral_flux_dn, "spectral_flux_dn_direct_sw");
      file.read(spectral_flux_up, "spectral_flux_up_sw");
      have_spectral_fluxes = true;
      have_band_fluxes = false;
    }
    else if (file.exist("band_flux_dn_direct_sw")) {
      file.read(spectral_flux_dn, "band_flux_dn_direct_sw");
      file.read(spectral_flux_up, "band_flux_up_sw");
      have_spectral_fluxes = true;
      have_band_fluxes = true;
      file.read(band_wavenumber1_, "band_wavenumber1_sw");
      file.read(band_wavenumber2_, "band_wavenumber2_sw");
    }

    if (have_spectral_fluxes) {
      spectral_flux_dn_.resize(ncol_new,nlev+1,spectral_flux_dn.size(3));
      spectral_flux_up_.resize(ncol_new,nlev+1,spectral_flux_up.size(3));
      int icol_new = 0;
      for (int icol = 0; icol < ncol; ++icol) {
	for (int isza = 0; isza < nsza; ++isza) {
	  spectral_flux_dn_(icol_new,__,__) = spectral_flux_dn(icol,index_sza(isza),__,__);
	  spectral_flux_up_(icol_new++,__,__) = spectral_flux_up(icol,index_sza(isza),__,__);
	}
      }

      effective_spectral_albedo_ = sum(spectral_flux_up_(__,end,__),0)
 	                         / sum(spectral_flux_dn_(__,end,__),0);

      if (have_band_fluxes && !band_mapping.empty()) {
	// Optionally map from narrow to wide bands
	int nband = maxval(band_mapping)+1;
	LOG << "  Mapping fluxes from " << band_mapping.size()
	    << " to " << nband << " bands\n";

	Vector new_band_wn1(nband), new_band_wn2(nband);
	Array3 new_spectral_flux_dn(ncol_new, nlev+1, nband);
	Array3 new_spectral_flux_up(ncol_new, nlev+1, nband);
	for (int jband = 0; jband < nband; ++jband) {
	  new_spectral_flux_dn(__,__,jband) = sum(spectral_flux_dn_(__,__,find(band_mapping==jband)),2);
	  new_spectral_flux_up(__,__,jband) = sum(spectral_flux_up_(__,__,find(band_mapping==jband)),2);
	}
	spectral_flux_dn_.clear();
	spectral_flux_dn_ = new_spectral_flux_dn;
	spectral_flux_up_.clear();
	spectral_flux_up_ = new_spectral_flux_up;
	effective_spectral_albedo_.clear();
	effective_spectral_albedo_ = sum(spectral_flux_up_(__,end,__),0)
	                           / sum(spectral_flux_dn_(__,end,__),0);
	
	for (int jband = 0; jband < nband; ++jband) {
	  new_band_wn1(jband) = minval(band_wavenumber1_(find(band_mapping==jband)));
	  new_band_wn2(jband) = maxval(band_wavenumber2_(find(band_mapping==jband)));
	}
	band_wavenumber1_.clear();
	band_wavenumber1_ = new_band_wn1;
	band_wavenumber2_.clear();
	band_wavenumber2_ = new_band_wn2;
      }
    }

    // Read high resolution fluxes at boundaries, if present
    if (file.exist("spectral_flux_dn_direct_surf_sw") 
	&& file.exist("spectral_flux_up_toa_sw")) {
      if (g_point.empty()) {
	WARNING << "Surface/TOA spectral fluxes ignored because g-point file not provided";
	ENDWARNING;
      }
      else {
	LOG << "  Mapping high-resolution boundary fluxes to g-points";
	int ng = maxval(g_point)+1;
	spectral_flux_dn_surf_.resize(ncol_new,ng);
	spectral_flux_up_toa_.resize(ncol_new,ng);
	
	int icol_new = 0;
	Vector spectral_flux_dn_surf, spectral_flux_up_toa;
	for (int icol = 0; icol < ncol; ++icol) {
	  LOG << ".";
	  for (int isza = 0; isza < nsza; ++isza) {
	    file.read(spectral_flux_up_toa, "spectral_flux_up_toa_sw", icol, index_sza(isza));
	    file.read(spectral_flux_dn_surf, "spectral_flux_dn_direct_surf_sw", icol, index_sza(isza));
	    for (int ig = 0; ig < ng; ++ig) {
	      intVector index = find(g_point == ig);
	      spectral_flux_dn_surf_(icol_new,ig) = sum(spectral_flux_dn_surf(index));
	      spectral_flux_up_toa_(icol_new,ig)  = sum(spectral_flux_up_toa(index));
	    }
	    ++icol_new;
	  }
	}
	LOG << "\n";
      }
    }

    ncol = ncol_new;
  }
  else {
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
    
    // Read high resolution fluxes at boundaries, if present
    if (file.exist("spectral_flux_dn_surf_lw") 
	&& file.exist("spectral_flux_up_toa_lw")) {
      if (g_point.empty()) {
	WARNING << "Surface/TOA spectral fluxes ignored because g-point file not provided";
	ENDWARNING;
      }
      else {
	LOG << "  Mapping high-resolution boundary fluxes to g-points";
	int ng = maxval(g_point)+1;
	spectral_flux_dn_surf_.resize(ncol,ng);
	spectral_flux_up_toa_.resize(ncol,ng);
	
	Vector spectral_flux_dn_surf, spectral_flux_up_toa;
	for (int icol = 0; icol < ncol; ++icol) {
	  LOG << ".";
	  file.read(spectral_flux_up_toa, "spectral_flux_up_toa_lw", icol);
	  file.read(spectral_flux_dn_surf, "spectral_flux_dn_surf_lw", icol);
	  for (int ig = 0; ig < ng; ++ig) {
	    intVector index = find(g_point == ig);
	    spectral_flux_dn_surf_(icol,ig) = sum(spectral_flux_dn_surf(index));
	    spectral_flux_up_toa_(icol,ig)  = sum(spectral_flux_up_toa(index));
	  }
	}
	LOG << "\n";
      }
    }

  }

  std::string molecules_str, molecule;
  file.read(molecules_str, DATA_FILE_GLOBAL_SCOPE, "constituent_id");
  std::stringstream molecules_s(molecules_str);
  LOG << "  Contains " << molecules_str << "\n";
  while (std::getline(molecules_s, molecule, ' ')) {
    std::size_t first_hyphen = molecule.find_first_of("-");
    if (first_hyphen == std::string::npos) {
      molecules_.push_back(molecule);
    }
    else {
      // If we have h2o-no-continuum then rename to h2o
      std::string new_molecule = molecule.substr(0,first_hyphen);
      LOG << "  Renaming " << molecule << " to " << new_molecule << "\n";
      molecules_.push_back(new_molecule);
    }
  }

  int nspec   = spectral_flux_up_.dimension(2);
  //  int nband= band_flux_up_.dimension(2);

  heating_rate_.resize(ncol,nlev);

  if (is_sw_) {
    // Heating rates neglect the upwelling contribution in the
    // shortwave, since it is affected by Rayleigh scattering that is
    // more expensive to compute in the iterative optimization
    for (int icol = 0; icol < ncol; ++icol) {
      heating_rate_single(pressure_hl_(icol,__), flux_dn_(icol,__), Vector(),
			  heating_rate_(icol,__));
    }

    if (have_spectral_fluxes) {
      spectral_heating_rate_.resize(ncol,nlev,nspec);
      for (int icol = 0; icol < ncol; ++icol) {
	heating_rate(pressure_hl_(icol,__), spectral_flux_dn_(icol,__,__), Matrix(),
		     spectral_heating_rate_(icol,__,__));
      }
    }
  }
  else {
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

void
LblFluxes::make_gas_mapping(const std::vector<std::string>& molecules)
{
  int ngas = molecules.size();
  gas_mapping.resize(ngas);
  gas_mapping = -1;
  for (int igas = 0; igas < ngas; ++igas) {
    for (int igas2 = 0; igas2 < molecules_.size(); ++igas2) {
      if (molecules_[igas2] == molecules[igas]) {
	gas_mapping[igas] = igas2;
      }
    }
  }
}

void
LblFluxes::mask_rayleigh_up(Real max_no_rayleigh_wavenumber)
{
  static bool done_message = false;
  intVector index = find(band_wavenumber2_ > max_no_rayleigh_wavenumber);
  if (!index.empty()) {
    effective_spectral_albedo_(index) = 0.0;
    spectral_flux_up_(__,__,index) = 0.0;
    flux_up_ = 0.0;
    if (!done_message) {
      LOG << "Ignoring upwelling for bands " << index << " because Rayleigh scattering not modelled\n";
      LOG << "  Effective albedo (up/direct_dn) in each bands is " << effective_spectral_albedo_ << "\n";
      done_message = true;
    }
  }
}

void
LblFluxes::subtract(const LblFluxes& source)
{
  flux_up_ -= source.flux_up_;
  flux_dn_ -= source.flux_dn_;
  spectral_flux_up_ -= source.spectral_flux_up_;
  spectral_flux_dn_ -= source.spectral_flux_dn_;
  heating_rate_ -= source.heating_rate_;
  spectral_heating_rate_ -= source.spectral_heating_rate_;
}

void
LblFluxes::calc_ckd_fluxes(const Array3D& optical_depth,
			   Array3D& flux_dn, Array3D& flux_up) const
{
  int nprof= optical_depth.size(0);
  int nlay = optical_depth.size(1);
  int ng   = optical_depth.size(2);

  flux_dn.resize(nprof,nlay+1,ng);
  flux_up.resize(nprof,nlay+1,ng);

  for (int iprof = 0; iprof < nprof; ++iprof) {
    if (is_sw_) {
      Real tsi_scaling = tsi_ / sum(solar_irradiance_);
      //      radiative_transfer_direct_sw(mu0_(iprof), tsi_scaling * solar_irradiance_, optical_depth(iprof,__,__),
      //				   flux_dn(iprof,__,__));
      Vector albedo_gpoint = effective_spectral_albedo_(iband_per_g);
      radiative_transfer_norayleigh_sw(mu0_(iprof), tsi_scaling * solar_irradiance_,
				     optical_depth(iprof,__,__), albedo_gpoint,
				     flux_dn(iprof,__,__), flux_up(iprof,__,__));
      //      WARNING << "LblFluxes::calc_ckd_fluxes setting upwelling solar fluxes to zero\n";
      //      ENDWARNING;
      //      flux_up(iprof,__,__) = 0.0;
    }
    else {
      radiative_transfer_lw(planck_hl_(iprof,__,__), optical_depth(iprof,__,__),
			    surf_emissivity_(iprof,iband_per_g), surf_planck_(iprof,__),
			    flux_dn(iprof,__,__), flux_up(iprof,__,__));
    }
  }
}
