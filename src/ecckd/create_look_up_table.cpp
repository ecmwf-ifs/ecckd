#include "read_spectrum.h"
#include "read_merged_spectrum.h"
#include "planck_function.h"
#include "OutputDataFile.h"
#include "ckd_model.h"
#include "constants.h"
#include "average_optical_depth.h"

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

  // Names of input and output file
  std::string input, output, ssi_file_name;

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  // Load g-point file
  if (!config.read(input, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  // The high-resolution solar spectral irradiance is used as a weight
  // when averaging gas optical depths
  Vector ssi;
  bool do_sw = false;
  if (config.read(ssi_file_name, "ssi")) {
    DataFile ssi_file(ssi_file_name);
    ssi_file.read(ssi, "solar_spectral_irradiance");
    do_sw = true;
  }

  intVector g_point;
  Vector band_wn1, band_wn2, solar_irradiance;
  intVector band_number;
  bool is_sw = false;
  std::string input_history, input_config;
  {
    LOG << "Reading " << input << "\n";
    DataFile input_file(input);
    if (!input_file.read(g_point, "g_point")) {
      ERROR << "\"g_point\" not found in \"" << input << "\"";
      THROW(PARAMETER_ERROR);
    }
    input_file.read(band_wn1, "wavenumber1_band");
    input_file.read(band_wn2, "wavenumber2_band");
    input_file.read(band_number, "band_number");
    // Solar irradiance in each g point copied from g-point definition
    // file to ckd definition file
    if (input_file.exist("solar_irradiance")) {
      input_file.read(solar_irradiance, "solar_irradiance");
      is_sw = true;
    }
    input_file.read(input_history, DATA_FILE_GLOBAL_SCOPE, "history");
    input_file.read(input_config, DATA_FILE_GLOBAL_SCOPE, "config");
  } // input_file is implicitly closed here

  int ng = maxval(g_point) + 1;

  // Check for and remove g points that refer to none of the spectrum
  std::vector<int> bad_g_points;
  for (int ig = 0; ig < ng; ++ig) {
    if (!any(g_point == ig)) {
      bad_g_points.push_back(ig);
    }
  }
  if (!bad_g_points.empty()) {
    LOG << "Removing " << bad_g_points.size() << " g point(s) that occupies none of the spectrum\n";
    int new_ng = ng-bad_g_points.size();
    intVector g_point_map(new_ng);
    int inewg = 0;
    int ibad = 0;
    for (int ig = 0; ig < ng; ++ig) {
      if (ibad >= bad_g_points.size() || ig != bad_g_points[ibad]) {
	g_point_map(inewg) = ig;
	++inewg;
      }
      else {
	++ibad;	
      }	
    }
    
    int nwav = g_point.size();
    intVector new_g_point(nwav);
    intVector new_band_number(new_ng);
    new_g_point = -1;
    for (int inewg = 0; inewg < new_ng; ++inewg) {
      LOG << "  Mapping new g point " << inewg
	  << " -> old g point " << g_point_map(inewg) << "\n";
      new_g_point.where(g_point == g_point_map(inewg)) = inewg;
      new_band_number(inewg) = g_point_map(inewg);
    }
    if (any(new_g_point < 0)) {
      ERROR << "Some unassigned spectral points after mapping";
      THROW(1);
    }
    g_point = new_g_point;
    band_number.clear();
    band_number = new_band_number;

    if (is_sw) {
      Vector new_solar_irradiance(new_ng);
      new_solar_irradiance = solar_irradiance(g_point_map);
      solar_irradiance.clear();
      solar_irradiance = new_solar_irradiance;
    }

    ng = new_ng;
  }

  std::vector<SingleGasData<false> > single_gas_data;

  int ngas = 0;
  int nlay = 0;
  int ncol = 0;

  std::string gas_str;

  Vector pressure_fl;
  Matrix temperature_fl; // Dimensioned (temperature,pressure)
  Vector wavenumber_cm_1, d_wavenumber_cm_1;
  Matrix optical_depth;
  Matrix planck_fl;
  int nwav;

  int temperature_stride = 1;
  config.read(temperature_stride, "temperature_stride");

  std::string averaging_method = "transmission";
  config.read(averaging_method, "averaging_method");

  // Loop over gases
  while (config.read(gas_str, "gases", ngas)) {
    std::string Gas = gas_str;
    std::transform(Gas.begin(), Gas.end(), Gas.begin(), ::toupper);
    LOG << "Creating look-up table for " << Gas << " (gas number " << ngas << ")\n";

    single_gas_data.push_back(SingleGasData<false>(gas_str));
    SingleGasData<false>& this_gas = single_gas_data[ngas];

    std::string conc_dependence_str;
    if (!config.read(conc_dependence_str, gas_str, "conc_dependence")) {
      ERROR << gas_str + ".conc_dependence not found in configuration";
      THROW(PARAMETER_ERROR);
    }

    if (conc_dependence_str == "none") {
      this_gas.conc_dependence = NONE;
    }
    else if (conc_dependence_str == "linear") {
      this_gas.conc_dependence = LINEAR;
    }
    else if (conc_dependence_str == "lut") {
      this_gas.conc_dependence = LUT;
    }
    else if (conc_dependence_str == "relative-linear") {
      this_gas.conc_dependence = RELATIVE_LINEAR;
    }
    else {
      ERROR << "conc_dependence \"" << conc_dependence_str << "\" not understood";
      THROW(PARAMETER_ERROR);
    }

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	std::string molecules;
	Real reference_surface_vmr = 1.0;
	ncol = 1;
	int icol = 0;
	int ngas_local = 0;

	Vector pressure_hl, temperature_hl;
	Matrix vmr_fl;

	while (icol < ncol) {
	  LOG << "  Reading temperature profile " << icol*temperature_stride << " for " << gas_str << "\n";
	  read_merged_spectrum(config, icol*temperature_stride, gas_str + ".",
			       pressure_hl, temperature_hl,
			       wavenumber_cm_1, d_wavenumber_cm_1, optical_depth,
			       molecules, vmr_fl, &ngas_local, &ncol);

	  ncol = (ncol+temperature_stride-1) / temperature_stride;

	  if (icol == 0) {
	    nlay = pressure_hl.size()-1;
	    pressure_fl = 0.5 * (pressure_hl(range(0,end-1)) + pressure_hl(range(1,end)));
	    this_gas.molar_abs.resize(ncol,nlay,ng);
	    this_gas.molar_abs = 0.0;
	    temperature_fl.resize(ncol,nlay);
	    this_gas.composite_molecules = molecules;
	    this_gas.composite_vmr = vmr_fl;
	  }

	  Vector t_x_p = temperature_hl * pressure_hl;
	  temperature_fl(icol,__) = 0.5 * (t_x_p(range(0,end-1)) + t_x_p(range(1,end))) / pressure_fl;

	  if (!do_sw) {
	    LOG << "  Computing Planck function\n";
	    nwav = wavenumber_cm_1.size();
	    planck_fl.resize(nlay,nwav);
	    planck_function(temperature_fl(icol,__), wavenumber_cm_1, d_wavenumber_cm_1,
			    planck_fl);
	    LOG << "  Planck-weighted averaging optical depths for each g point\n";

	    average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					     g_point, optical_depth, planck_fl, averaging_method,
					     this_gas.molar_abs(icol,__,__));
	  }
	  else {
	    LOG << "  Solar-spectrum-weighted averaging optical depths for each g point\n";
	    Matrix ssi_matrix = spread<0>(ssi,pressure_fl.size());
	    average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					     g_point, optical_depth, ssi_matrix, averaging_method,
					     this_gas.molar_abs(icol,__,__));
	  }

	  ++icol;
	}
      }
      break;
    case LINEAR:
    case RELATIVE_LINEAR:
      {
	std::string molecule;
	Real reference_surface_vmr;
	ncol = 1;
	int icol = 0;

	std::string file_name;
	if (!config.read(file_name, gas_str, "input")) {
	  ERROR << gas_str << ".input not found";
	  THROW(PARAMETER_ERROR);
	}

	Vector pressure_hl, temperature_hl, vmr_fl;

	if (this_gas.conc_dependence == RELATIVE_LINEAR) {
	  if (!config.read(this_gas.reference_vmr, gas_str, "reference_conc")) {
	    ERROR << gas_str << ".reference_conc must be provided if conc_dependence is relative-linear";
	    THROW(PARAMETER_ERROR);
	  }
	}

	while (icol < ncol) {
	  LOG << "  Reading temperature profile " << icol*temperature_stride << " from " << file_name << "\n";
	  read_spectrum(file_name, icol*temperature_stride, pressure_hl, temperature_hl,
			wavenumber_cm_1, d_wavenumber_cm_1, optical_depth,
			molecule, reference_surface_vmr, vmr_fl, &ncol);

	  ncol = (ncol+temperature_stride-1) / temperature_stride;

	  if (icol == 0) {
	    nlay = pressure_hl.size()-1;
	    pressure_fl = 0.5 * (pressure_hl(range(0,end-1)) + pressure_hl(range(1,end)));
	    this_gas.molar_abs.resize(ncol,nlay,ng);
	    this_gas.molar_abs = 0.0;
	    temperature_fl.resize(ncol,nlay);
	  }

	  Vector t_x_p = temperature_hl * pressure_hl;
	  temperature_fl(icol,__) = 0.5 * (t_x_p(range(0,end-1)) + t_x_p(range(1,end))) / pressure_fl;

	  if (!do_sw) {
	    LOG << "  Computing Planck function\n";
	    nwav = wavenumber_cm_1.size();
	    planck_fl.resize(nlay,nwav);
	    planck_function(temperature_fl(icol,__), wavenumber_cm_1, d_wavenumber_cm_1,
			    planck_fl);
	    LOG << "  Planck-weighted averaging optical depths for each g point\n";

	    average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					     g_point, optical_depth, planck_fl, averaging_method,
					     this_gas.molar_abs(icol,__,__));
	  }
	  else {
	    LOG << "  Solar-spectrum-weighted averaging optical depths for each g point\n";
	    Matrix ssi_matrix = spread<0>(ssi,pressure_fl.size());
	    average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					     g_point, optical_depth, ssi_matrix, averaging_method,
					     this_gas.molar_abs(icol,__,__));
	  }

	  ++icol;
	}

      }
      break;
    case LUT:
      {
	std::string molecule;
	Real reference_surface_vmr;
	ncol = 1;

	intVector nconc_vec = config.size(gas_str, "input");
	int nconc = nconc_vec[0];
	if (nconc <= 0) {
	  ERROR << gas_str << ".input not found";
	  THROW(PARAMETER_ERROR);
	}

	int iconc = 0;

	Vector pressure_hl, temperature_hl, vmr_fl;

	std::string file_name;
	while (config.read(file_name, gas_str, "input", iconc)) {
	  int icol = 0;
	  while (icol < ncol) {
	    LOG << "  Reading temperature profile " << icol*temperature_stride << " from " << file_name << "\n";
	    read_spectrum(file_name, icol*temperature_stride, pressure_hl, temperature_hl,
			  wavenumber_cm_1, d_wavenumber_cm_1, optical_depth,
			  molecule, reference_surface_vmr, vmr_fl, &ncol);

	    ncol = (ncol+temperature_stride-1) / temperature_stride;

	    if (iconc == 0 && icol == 0) {
	      nlay = pressure_hl.size()-1;
	      pressure_fl = 0.5 * (pressure_hl(range(0,end-1)) + pressure_hl(range(1,end)));
	      this_gas.molar_abs_conc.resize(nconc,ncol,nlay,ng);
	      this_gas.molar_abs_conc = 0.0;
	      this_gas.vmr.resize(nconc);
	      temperature_fl.resize(ncol,nlay);
	    }
	    this_gas.vmr(iconc) = reference_surface_vmr;

	    Vector t_x_p = temperature_hl * pressure_hl;
	    temperature_fl(icol,__) = 0.5 * (t_x_p(range(0,end-1)) + t_x_p(range(1,end))) / pressure_fl;

	    if (!do_sw) {
	      LOG << "  Computing Planck function\n";
	      nwav = wavenumber_cm_1.size();
	      planck_fl.resize(nlay,nwav);
	      planck_function(temperature_fl(icol,__), wavenumber_cm_1, d_wavenumber_cm_1,
			      planck_fl);
	      LOG << "  Planck-weighted averaging optical depths for each g point\n";

	      average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					       g_point, optical_depth, planck_fl, averaging_method,
					       this_gas.molar_abs_conc(iconc,icol,__,__));
	    }
	    else {
	      LOG << "  Solar-spectrum-weighted averaging optical depths for each g point\n";
	      Matrix ssi_matrix = spread<0>(ssi,pressure_fl.size());
	      average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					       g_point, optical_depth, ssi_matrix, averaging_method,
					       this_gas.molar_abs_conc(iconc,icol,__,__));
	    }

	    ++icol;
	  }
	  ++iconc;
	}
      }
    }


    ++ngas;
  }

  LOG << "Computing fraction of spectrum contributing to each g-point\n";

  // We store the fraction of the spectrum contributing to each
  // g-point in the variables gpoint_fraction, where the spectrum
  // intervals are bounded by wavenumber1 and wavenumber2.
  int dwav, startwav, endwav;

  if (!do_sw) {
    dwav = 10;
    startwav = 0;
    endwav = 3260;
  }
  else {
    dwav = 50;
    startwav = 250;
    endwav = 50000;
  }
  Vector wavenumber1 = static_cast<Real>(dwav) * range(startwav/dwav,endwav/dwav-1);
  Vector wavenumber2 = static_cast<Real>(dwav) * range(startwav/dwav+1,endwav/dwav);
  nwav = wavenumber1.size();

  Matrix gpoint_fraction(ng, nwav);
#pragma omp parallel for schedule (dynamic)
  for (int ig = 0; ig < ng; ++ig) {
    Real wav_per_gpoint = sum(d_wavenumber_cm_1.soft_link()(find(g_point == ig)));
    for (int iwav = 0; iwav < nwav; ++iwav) {
      gpoint_fraction(ig, iwav)
	= sum(d_wavenumber_cm_1.soft_link()(find(g_point == ig
				     && wavenumber_cm_1 >  wavenumber1(iwav)
				     && wavenumber_cm_1 <= wavenumber2(iwav))))
	/ wav_per_gpoint;
    }
  }

  LOG << "Writing " << output << "\n";

  std::string config_str;
  config.read(config_str);

  if (is_sw) {
    CkdModel<false> ckd_model(single_gas_data, solar_irradiance,
			      pressure_fl, temperature_fl,
			      wavenumber1, wavenumber2, gpoint_fraction,
			      band_wn1, band_wn2, band_number,
			      input_history, input_config);
    ckd_model.write(output, argc, argv, config_str);
  }
  else {

    LOG << "Generating Planck-function look-up table\n";

    Vector temperature_lut = range(120, 350);
    int nlut = temperature_lut.size();
    Matrix planck_lut(nlut,ng);
    
    for (int ig = 0; ig < ng; ++ig) {
      intVector index = find(g_point == ig);
      Matrix tmp_planck(nlut,index.size());
      planck_function(temperature_lut, wavenumber_cm_1(index), d_wavenumber_cm_1(index),
		      tmp_planck);
      planck_lut(__,ig) = sum(tmp_planck,1);
    }

    CkdModel<false> ckd_model(single_gas_data, temperature_lut, planck_lut,
			      pressure_fl, temperature_fl,
			      wavenumber1, wavenumber2, gpoint_fraction,
			      band_wn1, band_wn2, band_number,
			      input_history, input_config);
    ckd_model.write(output, argc, argv, config_str);
  }

  return 0;
}
