
#include "ckd_model.h"
#include "DataFile.h"
#include "file_manager.h"
#include "floating_point_exceptions.h"


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

  CkdModel<false> ckd_model(input);

  // Read location of g-points
  Vector wavenumber_cm_1;
  intVector g_point;
  int ng = ckd_model.ng();
  // First check if the g_points are stored in the raw CKD file - they
  // should be if they were modified from those in the find_g_points
  // program
  if (!ckd_model.read_g_points(wavenumber_cm_1, g_point)) {
    std::string gpoint_filename;
    if (!config.read(gpoint_filename, "gpointfile")) {
      ERROR << "gpointfile not provided";
      THROW(PARAMETER_ERROR);
    }
    DataFile gpointfile(gpoint_filename);
    gpointfile.read(wavenumber_cm_1, "wavenumber");
    gpointfile.read(g_point, "g_point");
    if (ng != maxval(g_point)+1) {
      ERROR << "Number of g-points in " << input << " does not match number in " << gpoint_filename;
      THROW(PARAMETER_ERROR);
    }
  }

  // Read spectral direct fluxes computed for reference atmosphere
  Matrix spectral_flux_dn;
  Matrix mole_fraction;
  Vector pressure_hl;
  Vector temperature_fl;
  std::string molecules_str;
  Real mu0;
  int imu0 = 0;
  {
    std::string lbl_filename;
    if (!config.read(lbl_filename, "lblfile")) {
      ERROR << "lblfile not provided";
      THROW(PARAMETER_ERROR);
    }
    LOG << "Reading " << lbl_filename << "\n";
    DataFile lblfile(lbl_filename);
    lblfile.read(mu0, "mu0", imu0);
    lblfile.read(molecules_str, DATA_FILE_GLOBAL_SCOPE, "constituent_id");
    lblfile.read(pressure_hl, "pressure_hl", imu0);
    Vector temperature_hl;
    lblfile.read(temperature_hl, "temperature_hl", imu0);
    temperature_fl = 0.5*(temperature_hl(range(0,end-1))+temperature_hl(range(1,end)));
    lblfile.read(mole_fraction, "mole_fraction_fl", imu0);
    lblfile.read(spectral_flux_dn, "spectral_flux_dn_direct_sw", imu0);
  }

  int nz = spectral_flux_dn.size(0)-1;
  int ngas = mole_fraction.size(0);

  LOG << "Computing optimal layer optical depths in each g point\n";
  Matrix od_best(nz, ng);
  for (int ig = 0; ig < ng; ++ig) {
    intVector index = find(g_point == ig);
    Real flux_top = sum(spectral_flux_dn(0, index));
    for (int iz = 0; iz < nz; ++iz) {
      Real flux_base = sum(spectral_flux_dn(iz+1,index));
      if (flux_base <= 0.0) {
	od_best(iz, ig) = -1.0;
      }
      else {
	od_best(iz, ig) = -mu0 * std::log(flux_base / flux_top);
      }
      flux_top = flux_base;
    }
  }

  

  LOG << "Running CKD model\n";
  spectral_flux_dn.clear();
  Matrix od_total(nz, ng);
  od_total = 0.0;
  std::stringstream molecules_s(molecules_str);
  for (int igas = -1; igas < ngas; ++igas) {
    std::string molecule;
    if (igas == -1) {
      molecule = "composite";
    }
    else {
      std::getline(molecules_s, molecule, ' ');
    }
    // calc_optical_depth requires arguments to have a column
    // dimension, but we are processing only one column so add a
    // singleton dimension to each argument with reshape
    int gas_index = ckd_model.get_gas_index(molecule);
    if (gas_index >= 0) {
      LOG << "  Gas " << igas << ": " << molecule << "\n";
      Array3D od;
      if (igas == -1) {
	od = ckd_model.calc_optical_depth(molecule,
					  pressure_hl.reshape(1,nz+1),
					  temperature_fl.reshape(1,nz));
      }
      else {
	od = ckd_model.calc_optical_depth(molecule,
					  pressure_hl.reshape(1,nz+1),
					  temperature_fl.reshape(1,nz),
					  mole_fraction[igas].reshape(1,nz));
      }
      od_total += od[0];
    }
    else {
      LOG << "  Gas " << igas << ": " << molecule << " not found\n";
    }
  }

  LOG << "Scaling coefficients in CKD look-up tables\n";
  Matrix scaling = od_best / od_total;
  scaling.where(od_best <= 0.0) = 1.0;
  Vector pressure_fl = 0.5*(pressure_hl(range(0,end-1))+pressure_hl(range(1,end)));
  ckd_model.scale_optical_depth(pressure_fl, scaling);

  std::string config_str;
  config.read(config_str);  

  ckd_model.write(output, argc, argv, config_str);

  return 0;
}
