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
  std::string input, output;

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  // Load g-point file
  if (!config.read(input, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  intVector g_point;

  {
    LOG << "Reading " << input << "\n";
    DataFile input_file(input);
    if (!input_file.read(g_point, "g_point")) {
      ERROR << "\"g_point\" not found in \"" << input << "\"";
      THROW(PARAMETER_ERROR);
    }
  } // input_file is implicitly closed here

  int ng = maxval(g_point) + 1;

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

	  LOG << "  Computing Planck function\n";
	  nwav = wavenumber_cm_1.size();
	  planck_fl.resize(nlay,nwav);
	  planck_function(temperature_fl(icol,__), wavenumber_cm_1, d_wavenumber_cm_1,
			  planck_fl);
	  LOG << "  Averaging optical depths for each g point\n";

	  average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					   g_point, optical_depth, planck_fl,
					   this_gas.molar_abs(icol,__,__));

	  ++icol;
	}
      }
      break;
    case LINEAR:
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

	  LOG << "  Computing Planck function\n";
	  nwav = wavenumber_cm_1.size();
	  planck_fl.resize(nlay,nwav);
	  planck_function(temperature_fl(icol,__), wavenumber_cm_1, d_wavenumber_cm_1,
			  planck_fl);
	  LOG << "  Averaging optical depths for each g point\n";

	  average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					   g_point, optical_depth, planck_fl,
					   this_gas.molar_abs(icol,__,__));

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

	    LOG << "  Computing Planck function\n";
	    nwav = wavenumber_cm_1.size();
	    planck_fl.resize(nlay,nwav);
	    planck_function(temperature_fl(icol,__), wavenumber_cm_1, d_wavenumber_cm_1,
			    planck_fl);
	    LOG << "  Averaging optical depths for each g point\n";

	    average_optical_depth_to_g_point(ng, reference_surface_vmr, pressure_fl, pressure_hl,
					     g_point, optical_depth, planck_fl,
					     this_gas.molar_abs_conc(iconc,icol,__,__));

	    ++icol;
	  }
	  ++iconc;
	}
      }
    }


    ++ngas;
  }


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

  LOG << "Computing fraction of spectrum contributing to each g-point\n";

  // We store the fraction of the spectrum contributing to each
  // g-point in the variables gpoint_fraction, where the spectrum
  // intervals are bounded by wavenumber1 and wavenumber2.
  int dwav = 10;
  int startwav = 0;
  int endwav = 3260;
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

  CkdModel<false> ckd_model(single_gas_data, temperature_lut, planck_lut,
			    pressure_fl, temperature_fl,
			    wavenumber1, wavenumber2, gpoint_fraction);
  ckd_model.write(output, argc, argv, config_str);

  /*

  OutputDataFile file(output);

  file.define_dimension("temperature", ncol);
  file.define_dimension("pressure", nlay);
  file.define_dimension("g_point", ng);
  file.define_dimension("temperature_planck", nlut);
  file.define_dimension("wavenumber", nwav);

  std::string molecule_list;
  for (int igas = 0; igas < ngas; ++igas) {
    if (igas == 0) {
      molecule_list = single_gas_data[igas].molecule;
    }
    else {
      molecule_list += " " + single_gas_data[igas].molecule;
    }
  }

  file.define_variable("n_gases", INT);
  file.write_long_name("Number of gases treated", "n_gases");
  file.write_comment("The gases are listed in the global attribute \"constituent_id\".", "n_gases");

  file.define_variable("temperature", FLOAT, "temperature", "pressure");
  file.write_long_name("Temperature", "temperature");
  file.write_units("K", "temperature");

  file.define_variable("pressure", FLOAT, "pressure");
  file.write_long_name("Pressure", "pressure");
  file.write_units("Pa", "pressure");

  file.define_variable("temperature_planck", FLOAT, "temperature_planck");
  file.write_long_name("Temperature for Planck function look-up table", "temperature_planck");
  file.write_units("K", "temperature_planck");
  file.define_variable("planck_function", FLOAT, "temperature_planck", "g_point");
  file.write_long_name("Planck function look-up table", "planck_function");
  file.write_units("W m-2", "planck_function");

  file.define_variable("wavenumber1", FLOAT, "wavenumber");
  file.write_long_name("Lower bound of spectral interval", "wavenumber1");
  file.write_units("cm-1", "wavenumber1");
  file.define_variable("wavenumber2", FLOAT, "wavenumber");
  file.write_long_name("Upper bound of spectral interval", "wavenumber2");
  file.write_units("cm-1", "wavenumber2");
  file.define_variable("gpoint_fraction", FLOAT, "g_point", "wavenumber");
  file.write_long_name("Fraction of spectrum contributing to each g-point",
		       "gpoint_fraction");

  std::string title = "Gas optics definition";
  file.write(title, "title");
  file.write("ecCKD gas optics tool", "source");

  file.write(molecule_list, "constituent_id");

  for (int igas = 0; igas < ngas; ++igas) {
    const SingleGasData<false>& this_gas = single_gas_data[igas];
    const std::string& Molecule = this_gas.Molecule;
    const std::string& molecule = this_gas.molecule;
    std::string varname = molecule + "_" + K_NAME;

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	file.define_variable(varname,
			     FLOAT, "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of background gases", varname);
	file.write_units("m2 mol-1", varname);
	file.write_comment("This is the absorption cross section of background gases per mole of dry air.",
			   varname);
	file.define_dimension(molecule + "_gas", this_gas.composite_vmr.size(0));
	file.define_variable(molecule + "_mole_fraction", FLOAT, molecule + "_gas", "pressure");
	file.write_long_name("Mole fractions of the gases that make up " + Molecule, 
			     molecule + "_mole_fraction");
	file.write_units("1", molecule + "_mole_fraction");
	file.write_comment("The gases that make up " + Molecule + " are listed in the global attribute \""
		   + molecule + "_constituent_id\".", molecule + "_mole_fraction");
	file.write(this_gas.composite_molecules, molecule + "_constituent_id");
      }
      break;
    case LINEAR:
      {
	file.define_variable(varname,
			     FLOAT, "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of " + Molecule,
			     varname);
	file.write_units("m2 mol-1", varname);
      }
      break;
    case LUT:
      {
	file.define_dimension(molecule + "_mole_fraction", this_gas.vmr.size());

	file.define_variable(molecule + "_mole_fraction", FLOAT, molecule + "_mole_fraction");
	file.write_long_name(Molecule + " mole fraction for look-up table", molecule + "_mole_fraction");
	file.write_units("1", molecule + "_mole_fraction");

	file.define_variable(varname,
			     FLOAT, molecule + "_mole_fraction", "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of " + Molecule,
			     varname);
	file.write_units("m2 mol-1", varname);
      }
    }
  }

  file.append_history(argc, argv);

  std::string config_str;
  config.read(config_str);  
  file.write(config_str, "config");

  // Write data

  file.end_define_mode();

  file.write(ngas, "n_gases");
  file.write(pressure_fl, "pressure");
  file.write(temperature_fl, "temperature");
  file.write(temperature_lut, "temperature_planck");
  file.write(wavenumber1, "wavenumber1");
  file.write(wavenumber2, "wavenumber2");
  file.write(gpoint_fraction, "gpoint_fraction");

  for (int igas = 0; igas < ngas; ++igas) {
    const SingleGasData<false>& this_gas = single_gas_data[igas];
    const std::string& molecule = this_gas.molecule;

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	file.write(this_gas.composite_vmr, molecule + "_mole_fraction");
      }
      // No break here: fall through to next one
    case LINEAR:
      {
	file.write(this_gas.molar_abs, molecule + "_" + K_NAME);
      }
      break;
    case LUT:
      {
	file.write(this_gas.vmr, molecule + "_mole_fraction");
	for (int iconc = 0; iconc < this_gas.vmr.size(); ++iconc) {
	  file.write(this_gas.molar_abs_conc(iconc,__,__,__),
		     molecule + "_" + K_NAME, iconc);
	}
      }
    }
  }

  file.write(planck_lut, "planck_function");

  file.close();
  */
}
