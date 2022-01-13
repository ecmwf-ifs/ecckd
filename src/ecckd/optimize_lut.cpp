// optimize_lut.cpp - Optimize a gas-optics look-up table
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

#include "radiative_transfer_lw.h"
#include "ckd_model.h"
#include "lbl_fluxes.h"
#include "solve_adept.h"
#include "floating_point_exceptions.h"
#include "DataFile.h"
#include "file_manager.h"

// Main program 
int
main(int argc, const char* argv[])
{

  // CONFIGURATION

  // Read configuration information from command-line and first file
  // on command-line
  DataFile config(argc, argv);

  std::string log_level;
  if (config.read(log_level, "log_level")) {
    set_log_level(log_level);
  }

  // Modify search path
  std::string mod_path;
  if (config.read(mod_path, "prepend_path")) {
    prepend_search_directory(mod_path);
  }
  if (config.read(mod_path, "append_path")) {
    append_search_directory(mod_path);
  }

  // Allow debugging
  //  set_trace_exceptions(true);

  enable_floating_point_exceptions();

  // Name of input and output coefficient files
  std::string input, output;

  if (!config.read(input, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  // Gases to optimize (others kept constant)
  std::string gas;
  std::vector<std::string> gas_list;
  int igas = 0;

  LOG << "Optimizing coefficients of:";
  while (config.read(gas, "gases", igas)) {
    LOG << " " << gas;
    gas_list.push_back(gas);
    ++igas;
  }
  if (igas == 0) {
    LOG << " ALL GASES\n";
  }
  else {
    LOG << "\n";
  }

  Real flux_weight = 0.02;
  Real flux_profile_weight = 0.0;
  Real broadband_weight = 0.5;
  Real prior_error = 1.0;
  Real rayleigh_prior_error = 0.0;
  Real spectral_boundary_weight = 0.0;

  Real temperature_corr = 0.5;
  Real pressure_corr = 0.5;
  Real conc_corr = 0.5;
  Real convergence_criterion = 0.02;

  // Wavenumber above which is it not safe to neglect Rayleigh
  // scattering in calculating upwelling fluxes
  Real max_no_rayleigh_wavenumber = 10000.0;

  std::string model_id;

  // We may wish to merge some of the narrow bands in the LBL flux
  // files
  intVector band_mapping;

  config.read(flux_weight, "flux_weight");
  config.read(flux_profile_weight, "flux_profile_weight");
  config.read(broadband_weight, "broadband_weight");
  config.read(spectral_boundary_weight, "spectral_boundary_weight");
  config.read(prior_error, "prior_error");
  if (config.read(rayleigh_prior_error, "rayleigh_prior_error")) {
    if (rayleigh_prior_error > 0.0) {
      LOG << "Optimizing Rayleigh scattering coefficients with prior error of " << rayleigh_prior_error << "\n";
      gas_list.push_back("rayleigh");
    }
  }
  config.read(temperature_corr, "temperature_corr");
  config.read(pressure_corr, "pressure_corr");
  config.read(conc_corr, "conc_corr");
  config.read(convergence_criterion, "convergence_criterion");
  config.read(model_id, "model_id");
  config.read(max_no_rayleigh_wavenumber, "max_no_rayleigh_wavenumber");

  int max_iterations = 3000;
  config.read(max_iterations, "max_iterations");

  Real negative_od_penalty = 1.0e4;
  config.read(negative_od_penalty, "negative_od_penalty");

  if (config.exist("band_mapping")) {
    config.read(band_mapping, "band_mapping");
  }

  Stack stack;

  CkdModel<true> ckd_model(input, gas_list);
  ckd_model.set_model_id(model_id);

  // Optionally read location of g-points, since some training data
  // might contain high resolution boundary fluxes that need to be
  // mapped to g-points
  Vector wavenumber_cm_1;
  intVector g_point;
  int ng = ckd_model.ng();
  // First check if the g_points are stored in the raw CKD file - they
  // should be if they were modified from those in the find_g_points
  // program
  if (!ckd_model.read_g_points(wavenumber_cm_1, g_point)) {
    std::string gpoint_filename;
    if (config.read(gpoint_filename, "gpointfile")) {
      DataFile gpointfile(gpoint_filename);
      gpointfile.read(g_point, "g_point");
      if (ng != maxval(g_point)+1) {
	ERROR << "Number of g-points in " << input << " does not match number in " << gpoint_filename;
	THROW(PARAMETER_ERROR);
      }
    }
  }

  ckd_model.cap_relative_linear_coeffts(0.8);

  ckd_model.create_error_covariances(prior_error, pressure_corr, temperature_corr, conc_corr,
				     rayleigh_prior_error);

  // Optional: compute radiative transfer of one set of profiles
  // relative to another, useful to get the forcing of minor gases
  // correct
  std::string relative_to_file;
  bool have_relative_to_fluxes = false;
  LblFluxes relative_to_fluxes;

  Array3D* relative_ckd_flux_dn = 0;
  Array3D* relative_ckd_flux_up = 0;
  Array3D  rel_ckd_flux_dn;
  Array3D  rel_ckd_flux_up;

  if (config.read(relative_to_file, "relative_to")) {
    have_relative_to_fluxes = true;
    LOG << "Errors evaluated relative to the following file:\n";
    relative_to_fluxes.read(relative_to_file, band_mapping);

    relative_to_fluxes.make_gas_mapping(ckd_model.molecules);
    if (!ckd_model.is_sw()) {
      relative_to_fluxes.planck_hl_   = ckd_model.calc_planck_function(relative_to_fluxes.temperature_hl_);
      relative_to_fluxes.surf_planck_ = ckd_model.calc_planck_function(relative_to_fluxes.temperature_hl_(__,end));
    }
    else {
      relative_to_fluxes.solar_irradiance_ = ckd_model.solar_irradiance();
    }

    if (relative_to_fluxes.have_band_fluxes) {
      relative_to_fluxes.iband_per_g = ckd_model.iband_per_g(relative_to_fluxes.band_wavenumber1_,
							     relative_to_fluxes.band_wavenumber2_);
    }

    ADEPT_ACTIVE_STACK->new_recording();
    aArray3D aod;
    calc_total_optical_depth(ckd_model, relative_to_fluxes, aod, true);
    Array3D od = value(aod);
    relative_to_fluxes.calc_ckd_fluxes(od, rel_ckd_flux_dn, rel_ckd_flux_up);
    relative_ckd_flux_dn = &rel_ckd_flux_dn;
    relative_ckd_flux_up = &rel_ckd_flux_up;
  }

  // LBL fluxes used for training
  std::string training_file;

  std::vector<LblFluxes> training_data;

  int itrain = 0;
  while (config.read(training_file, "training_input", itrain)) {

    LblFluxes fluxes(training_file, band_mapping, g_point);
    if (have_relative_to_fluxes) {
      LOG << "  Subtracting reference fluxes\n";
      fluxes.subtract(relative_to_fluxes);
    }

    fluxes.make_gas_mapping(ckd_model.molecules);

    if (band_mapping.empty()) {
      if (fluxes.nspec() != ckd_model.ng()) {
	ERROR << "band_mapping not provided, so number of g-points must match between LBL and CKD models";
	THROW(PARAMETER_ERROR);
      }
    }

    if (!ckd_model.is_sw()) {
      fluxes.planck_hl_   = ckd_model.calc_planck_function(fluxes.temperature_hl_);
      fluxes.surf_planck_ = ckd_model.calc_planck_function(fluxes.temperature_hl_(__,end));
    }
    else {
      fluxes.solar_irradiance_ = ckd_model.solar_irradiance();
    }

    if (fluxes.have_band_fluxes) {
      fluxes.iband_per_g = ckd_model.iband_per_g(fluxes.band_wavenumber1_,
						 fluxes.band_wavenumber2_);
    }

    if (ckd_model.is_sw()) {
      // Remove upwelling fluxes affected by Rayleigh scattering since
      // they will not be adequately modelled by the fast scheme used
      // for the optimization
      fluxes.mask_rayleigh_up(max_no_rayleigh_wavenumber);
    }

    training_data.push_back(fluxes);

    ++itrain;
  }

  if (itrain == 0) {
    ERROR << "\"training_input\" not specified";
    THROW(PARAMETER_ERROR);
  }

  adept::MinimizerStatus status = solve_adept(ckd_model, training_data,
			   flux_weight, flux_profile_weight, broadband_weight,
			   spectral_boundary_weight, prior_error,
			   max_iterations, convergence_criterion,
			   negative_od_penalty,
			   relative_ckd_flux_dn, relative_ckd_flux_up);

  LOG << "Convergence status: " << adept::minimizer_status_string(status) << "\n";

  std::string config_str;
  config.read(config_str);  

  //  ckd_model.cap_relative_linear_coeffts(1.0);
  ckd_model.write(output, argc, argv, config_str);

  return 0;
}
