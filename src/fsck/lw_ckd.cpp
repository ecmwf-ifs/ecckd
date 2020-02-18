#include "DataFile.h"
#include "OutputDataFile.h"
#include "ckd_model.h"
#include "radiative_transfer_lw.h"

using namespace adept;

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

  std::string ckd_file, input_file, output_file;

  if (!config.read(ckd_file, "ckd_model")) {
    ERROR << "\"ckd_model\" not specified";
    THROW(PARAMETER_ERROR);
  }
  if (!config.read(input_file, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(output_file, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  // Gases to use
  std::string gas;
  std::vector<std::string> gas_list;
  int igas = 0;
  while (config.read(gas, "gases", igas)) {
    gas_list.push_back(gas);
    ++igas;
  }

  Real co2_scaling = -1.0;
  Real ch4_scaling = -1.0;
  Real n2o_scaling = -1.0;
  Real cfc11_scaling = -1.0;
  Real cfc12_scaling = -1.0;
  config.read(co2_scaling, "co2_scaling");
  config.read(ch4_scaling, "ch4_scaling");
  config.read(n2o_scaling, "n2o_scaling");
  config.read(cfc11_scaling, "cfc11_scaling");
  config.read(cfc12_scaling, "cfc12_scaling");

  CkdModel<false> ckd_model(ckd_file);

  DataFile input(input_file);

  Matrix temperature_hl, pressure_hl;

  input.read(temperature_hl, "temperature_hl");
  input.read(pressure_hl,    "pressure_hl");

  Matrix temperature_fl, p_x_t;
  p_x_t = temperature_hl * pressure_hl;
  temperature_fl = (p_x_t(__,range(0,end-1)) + p_x_t(__,range(1,end)))
    / (pressure_hl(__,range(0,end-1)) + pressure_hl(__,range(1,end)));

  Vector temperature_surf;
  temperature_surf = temperature_hl(__,end);

  Array3 planck_hl   = ckd_model.calc_planck_function(temperature_hl);
  Matrix planck_surf = ckd_model.calc_planck_function(temperature_surf);

  int ncol = temperature_hl.dimension(0);
  int np   = temperature_hl.dimension(1)-1;
  int ng   = ckd_model.ng();

  LOG << "Writing " << output_file << "\n";

  OutputDataFile file(output_file);

  bool write_od = true;

  // Define dimensions
  file.define_dimension("column", ncol);
  if (write_od) {
    file.define_dimension("level", np);
  }
  file.define_dimension("half_level", np+1);
  if (write_od) {
    file.define_dimension("g_point", ng);
  }

  // Define variables
  file.define_variable("pressure_hl", FLOAT, "column", "half_level");
  file.write_long_name("Pressure", "pressure_hl");
  file.write_units("Pa", "pressure_hl");

  if (write_od) {
    file.define_variable("optical_depth", FLOAT, "column", "level", "g_point");
    file.write_long_name("Optical depth", "optical_depth");

    for (int igas = 0; igas < ckd_model.molecules.size(); ++igas) {
      std::string& molecule = ckd_model.molecules[igas];
      file.define_variable(molecule + "_optical_depth", FLOAT, "column", "level", "g_point");
      file.write_long_name(molecule + " optical depth", molecule + "_optical_depth");
    }

    file.define_variable("planck_hl", FLOAT, "column", "half_level", "g_point");
    file.write_long_name("Planck function", "planck_hl");
    file.write_units("W m-2", "planck_hl");

    file.define_variable("planck_surf", FLOAT, "column", "g_point");
    file.write_long_name("Planck function at surface", "planck_surf");
    file.write_units("W m-2", "planck_surf");

    file.define_variable("spectral_flux_up_lw", FLOAT, "column", "half_level", "g_point");
    file.write_long_name("Spectral upwelling longwave flux", "spectral_flux_up_lw");
    file.write_units("W m-2", "spectral_flux_up_lw");

    file.define_variable("spectral_flux_dn_lw", FLOAT, "column", "half_level", "g_point");
    file.write_long_name("Spectral downwelling longwave flux", "spectral_flux_dn_lw");
    file.write_units("W m-2", "spectral_flux_dn_lw");
  }

  file.define_variable("flux_up_lw", FLOAT, "column", "half_level");
  file.write_long_name("Upwelling longwave flux", "flux_up_lw");
  file.write_units("W m-2", "flux_up_lw");

  file.define_variable("flux_dn_lw", FLOAT, "column", "half_level");
  file.write_long_name("Downwelling longwave flux", "flux_dn_lw");
  file.write_units("W m-2", "flux_dn_lw");
 
  // Write data
  file.end_define_mode();

  file.write(pressure_hl, "pressure_hl");

  Array3 od(ncol,np,ng), od_gas(ncol,np,ng);

  od = 0.0;

  for (int igas = 0; igas < ckd_model.molecules.size(); ++igas) {
    std::string& molecule = ckd_model.molecules[igas];

    // If gas list provided, check molecule is in it
    if (!gas_list.empty()) {
      if (std::find(gas_list.begin(), gas_list.end(), molecule) == gas_list.end()) {
	LOG << "  Skipping " << molecule << "\n";
	continue;
      }
    }

    Matrix vmr;
    std::string var_name = molecule + "_vmr";
    if (!input.exist(var_name)) {
      WARNING << "  " << var_name << " not found";
      ENDWARNING;
    }
    else {
      LOG << "  Computing optical depth of " << molecule;
      input.read(vmr, var_name);

      if (co2_scaling >= 0.0 && molecule == "co2") {
	vmr *= co2_scaling;
	LOG << " from concentration scaled by " << co2_scaling;
      }
      else if (ch4_scaling >= 0.0 && molecule == "ch4") {
	vmr *= ch4_scaling;
	LOG << " from concentration scaled by " << ch4_scaling;
      }
      else if (n2o_scaling >= 0.0 && molecule == "n2o") {
	vmr *= n2o_scaling;
	LOG << " from concentration scaled by " << n2o_scaling;
      }
      else if (cfc11_scaling >= 0.0 && molecule == "cfc11") {
	vmr *= cfc11_scaling;
	LOG << " from concentration scaled by " << cfc11_scaling;
      }
      else if (cfc12_scaling >= 0.0 && molecule == "cfc12") {
	vmr *= cfc12_scaling;
	LOG << " from concentration scaled by " << cfc12_scaling;
      }
      LOG << "\n";

      od_gas = ckd_model.calc_optical_depth(igas, pressure_hl, temperature_fl, vmr);
      od += od_gas;
      if (write_od) {
	file.write(od_gas, molecule + "_optical_depth");
      }
    }
  }

  Array3D flux_up(ncol,np+1,ng), flux_dn(ncol,np+1,ng);

  Vector surf_emissivity(ng);
  surf_emissivity = 1.0;

  
  for (int icol = 0; icol < ncol; ++icol) {
    radiative_transfer_lw(planck_hl(icol,__,__), od(icol,__,__), surf_emissivity,
			  planck_surf(icol,__), flux_dn(icol,__,__), flux_up(icol,__,__));
  }
 
  if (write_od) {
    file.write(od, "optical_depth");
    file.write(planck_hl, "planck_hl");
    file.write(planck_surf, "planck_surf");
  }

  file.write(flux_up, "spectral_flux_up_lw");
  file.write(flux_dn, "spectral_flux_dn_lw");

  file.write(sum(flux_up,2), "flux_up_lw");
  file.write(sum(flux_dn,2), "flux_dn_lw");

  file.close();

}
