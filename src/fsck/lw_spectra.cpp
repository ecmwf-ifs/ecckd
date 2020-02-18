#include <string>
#include <algorithm>

#include "read_merged_spectrum.h"
#include "planck_function.h"
#include "radiative_transfer_lw.h"
#include "average_optical_depth.h"
#include "OutputDataFile.h"

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

  // Name of output file
  std::string output;

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  std::string input;
  config.read(input, "input");

  DataFile gpoints;
  bool have_gpoints = false;
  Vector g_point;
  int ng = -1;
  if (config.exist("gpoints")) {
    std::string gpoint_file;
    config.read(gpoint_file, "gpoints");
    gpoints.open(gpoint_file);
    gpoints.read(g_point, "g_point");
    ng = maxval(g_point) + 1;
    have_gpoints = true;
  }

  // Index of profiles to extract from files
  int iprofile = -1;
  bool do_one_profile = config.read(iprofile, "iprofile");

  Matrix optical_depth;
    
  // Wavenumber, cm-1
  Vector wavenumber_cm_1;
  // Wavenumber spacing, cm-1
  Vector d_wavenumber_cm_1;
  // Half-level pressure (Pa) and temperature (K)
  Vector pressure_hl, temperature_hl;
  // Volume mixing ratio of each gas (mol mol-1)
  Matrix vmr_fl;
  std::string molecules;

  int ngas;
  int ncol = 10000;

  int icol = 0;

  if (do_one_profile) {
    icol = iprofile;
  }

  bool is_first_profile = true;

  OutputDataFile file(output);

  while (icol < ncol) {
    LOG << "Profile " << icol << "\n";
    read_merged_spectrum(config, icol, "",
			 pressure_hl, temperature_hl,
			 wavenumber_cm_1, d_wavenumber_cm_1,
			 optical_depth, molecules, vmr_fl, &ngas, &ncol);
    int nlay = optical_depth.size(0);
    int nwav = optical_depth.size(1);

    if (is_first_profile) {

      is_first_profile = false;

      int nprofile = ncol;
      if (do_one_profile) {
	iprofile = 0;
	nprofile = 1;
      }
      else {
	iprofile = icol;
      }

      // Define dimensions

      // If "column" is an unlimited dimension then it helps to use
      // the NCO tool "ncrcat" to concatenate multiple files.
      //      file.define_dimension("column", nprofile);
      file.define_dimension("column", 0);

      file.define_dimension("level", nlay);
      file.define_dimension("half_level", nlay+1);
      std::string spec_name;
      if (!have_gpoints) {
	spec_name = "wavenumber";
	file.define_dimension(spec_name, nwav);
      }
      else {
	spec_name = "g_point";
	file.define_dimension(spec_name, ng);
      }
      file.define_dimension("gas", ngas);

      // Define variables
      file.define_variable("pressure_hl", FLOAT, "column", "half_level");
      file.write_long_name("Pressure at half levels", "pressure_hl");
      file.write_units("Pa", "pressure_hl");

      file.define_variable("temperature_hl", FLOAT, "column", "half_level");
      file.write_long_name("Temperature at half levels", "temperature_hl");
      file.write_units("K", "temperature_hl");

      if (!have_gpoints) {
	file.define_variable("wavenumber", DOUBLE, "wavenumber");
	file.deflate_variable("wavenumber");
	file.write_long_name("Wavenumber", "wavenumber");
	file.write_units("cm-1", "wavenumber");
      }

      file.define_variable("vmr_fl", FLOAT, "column", "gas", "level");
      file.write_long_name("Volume mixing ratio", "vmr_fl");
      file.write_units("mol mol-1", "vmr_fl");
      file.write_comment("The gases are listed in the global attribute \"molecules\".", "vmr_fl");

      file.define_variable("flux_dn_lw", FLOAT, "column", "half_level");
      file.write_long_name("Upwelling longwave flux", "flux_dn_lw");
      file.write_units("W m-2", "flux_dn_lw");

      file.define_variable("flux_up_lw", FLOAT, "column", "half_level");
      file.write_long_name("Upwelling longwave flux", "flux_up_lw");
      file.write_units("W m-2", "flux_up_lw");

      file.define_variable("optical_depth", FLOAT, "column", "level", spec_name);
      if (!have_gpoints) {
	file.deflate_variable("optical_depth");
	file.set_chunking({1,nlay,nwav}, "optical_depth");
      }
      file.write_long_name("Layer optical depth", "optical_depth");

      file.define_variable("spectral_flux_dn_lw", FLOAT, "column", "half_level", spec_name);
      file.write_long_name("Downwelling longwave spectral flux", "spectral_flux_dn_lw");
      file.write_units("W m-2", "spectral_flux_dn_lw");

      file.define_variable("spectral_flux_up_lw", FLOAT, "column", "half_level", spec_name);
      file.write_long_name("Upwelling longwave spectral flux", "spectral_flux_up_lw");
      file.write_units("W m-2", "spectral_flux_up_lw");

      file.append_history(argc, argv);
      // Replace commas with spaces in "molecules"
      std::replace(molecules.begin(), molecules.end(), ',', ' ');
      file.write(molecules, "molecules");
      std::string config_str;
      config.read(config_str);  
      file.write(config_str, "config");

      file.end_define_mode();

    }

    file.write(pressure_hl, "pressure_hl", iprofile);
    file.write(temperature_hl, "temperature_hl", iprofile);
    
    file.write(vmr_fl, "vmr_fl", iprofile);

    LOG << "  Computing Planck function\n";
    Matrix planck_hl;
    Vector surf_planck, surf_emissivity;
    
    planck_hl.resize(nlay+1,nwav);
    planck_function(temperature_hl, wavenumber_cm_1, d_wavenumber_cm_1,
		    planck_hl);
    surf_planck.resize(nwav);
    planck_function(temperature_hl(end), wavenumber_cm_1, d_wavenumber_cm_1,
		    surf_planck);
    surf_emissivity.resize(nwav);
    surf_emissivity = 1.0;
    
    Matrix flux_dn(nlay+1,nwav), flux_up(nlay+1,nwav);
    
    LOG << "  Performing longwave radiative transfer\n";
    radiative_transfer_lw(planck_hl, optical_depth, surf_emissivity,
			  surf_planck, flux_dn, flux_up);
    
    planck_hl.clear();
    surf_planck.clear();
    surf_emissivity.clear();
    
    file.write(sum(flux_dn,1), "flux_dn_lw", iprofile);
    file.write(sum(flux_up,1), "flux_up_lw", iprofile);
      
    if (!have_gpoints) {
      file.write(optical_depth, "optical_depth", iprofile);
      file.write(flux_dn, "spectral_flux_dn_lw", iprofile);
      file.write(flux_up, "spectral_flux_up_lw", iprofile);
    }
    else {
      Matrix spectral_flux_dn(nlay+1,ng), spectral_flux_up(nlay+1,ng);
      Matrix spectral_od(nlay,ng);
      Vector pressure_fl = 0.5 * (pressure_hl(range(0,end-1)) + pressure_hl(range(1,end)));
      Vector t_x_p = temperature_hl * pressure_hl;
      Vector temperature_fl = 0.5 * (t_x_p(range(0,end-1)) + t_x_p(range(1,end))) / pressure_fl;
      Matrix planck_fl(nlay,nwav);

      planck_function(temperature_fl, wavenumber_cm_1, d_wavenumber_cm_1,
		      planck_fl);
      average_optical_depth_to_g_point(ng, 0.0, pressure_fl, pressure_hl,
				       g_point, optical_depth, planck_fl,
				       spectral_od);
      for (int ig = 0; ig < ng; ++ig) {
	intVector index = find(g_point == ig);

	spectral_flux_dn(__,ig) = sum(flux_dn(__,index),1);
	spectral_flux_up(__,ig) = sum(flux_up(__,index),1);
      }
      file.write(spectral_od, "optical_depth", iprofile);
      file.write(spectral_flux_dn, "spectral_flux_dn_lw", iprofile);
      file.write(spectral_flux_up, "spectral_flux_up_lw", iprofile);
    }

    if (do_one_profile) {
      break;
    }
    ++icol;
    ++iprofile;
  }

  file.close();

}
