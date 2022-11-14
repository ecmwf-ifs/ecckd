// ckd_model.cpp - Class for storing a CKD model
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

#include <algorithm>
#include <sstream>
#include <cmath>
#include <iostream>

#include "ckd_model.h"
#include "DataFile.h"
#include "OutputDataFile.h"
#include "constants.h"
#include "adept_scalar.h"
#include "write_standard_attributes.h"
#include "rayleigh_scattering.h"

using namespace adept;

template<bool IsActive>
void
CkdModel<IsActive>::read(const std::string& file_name, 
			 const std::vector<std::string>& active_gas_list)
{
  clear();
  LOG << "Reading CKD definition file " << file_name << "\n";
  DataFile file(file_name);
  file.throw_exceptions(true);
  if (file.exist("solar_irradiance")) {
    file.read(solar_irradiance_, "solar_irradiance");
    if (file.exist("solar_spectral_irradiance")) {
      file.read(ssi_, "solar_spectral_irradiance");
    }
  }
  else {
    file.read(temperature_planck_, "temperature_planck");
    file.read(planck_function_, "planck_function");
  }
  file.read(temperature_, "temperature");

  {
    Vector pressure;
    file.read(pressure, "pressure");
    log_pressure_ = log(pressure);
  }

  file.read(wavenumber1_, "wavenumber1");
  file.read(wavenumber2_, "wavenumber2");
  file.read(gpoint_fraction_, "gpoint_fraction");
  file.read(wavenumber1_band_, "wavenumber1_band");
  file.read(wavenumber2_band_, "wavenumber2_band");
  file.read(band_number_, "band_number");

  // Optional input
  if (file.exist("g_point")) {
    file.read(wavenumber_hr_, "wavenumber_hr");
    file.read(g_point_, "g_point");
  }

  np_   = log_pressure_.size();
  nt_   = temperature_.dimension(0);
  ng_   = gpoint_fraction_.dimension(0);
  nwav_ = gpoint_fraction_.dimension(1);

  std::string molecules_str;
  file.read(molecules_str, DATA_FILE_GLOBAL_SCOPE, "constituent_id");

  file.read(history_, DATA_FILE_GLOBAL_SCOPE, "history");
  file.read(summary_, DATA_FILE_GLOBAL_SCOPE, "summary");
  file.read(config_,  DATA_FILE_GLOBAL_SCOPE, "config");
  file.read(model_id_,DATA_FILE_GLOBAL_SCOPE, "model_id");

  int n_gases;
  file.read(n_gases, "n_gases");

  std::string molecule;
  int igas = 0;

  // Count number of state variables
  int nx = 0;

  {
    std::stringstream molecules_s(molecules_str);
    while (std::getline(molecules_s, molecule, ' ')) {
      intVector ndims = file.size(molecule + "_" + K_NAME);
      if (is_active_(std::string(molecule), active_gas_list)) {
	nx += product(ndims);
      }
    }
  }

  if (is_sw() && is_active_("rayleigh", active_gas_list)) {
    nx += ng_;
  }

  if (IsActive) {
    x.resize(nx);
    LOG << "  State variable size: " << nx << "\n";
  }
  int ix = 0; // Current index to state vector x

  bool is_first = true;
  bool have_min_max = false;

  // Read in molar absorption coefficients
  std::stringstream molecules_s(molecules_str);
  while (std::getline(molecules_s, molecule, ' ')) {
    single_gas_data_.push_back(SingleGasData<IsActive>(molecule));
    molecules.push_back(molecule);
    SingleGasData<IsActive>& this_gas = single_gas_data_[igas];

    LOG << "  Reading absorption properties of " << this_gas.Molecule << "\n";

    // First check if we have min/max values
    if (is_first) {
      if (file.exist(molecule + "_" + K_NAME + "_min")) {
	have_min_max = true;
	if (IsActive) {
	  x_min.resize(nx);
	  x_max.resize(nx);
	}
      }
      is_first = false;
    }

    if (file.exist(molecule + "_mole_fraction")
	&& file.size(molecule + "_mole_fraction").size() == 1) {
      this_gas.conc_dependence = LUT;
      LOG << "    Concentration dependence: look-up table\n";
      file.read(this_gas.vmr, molecule + "_mole_fraction");
      int nconc = this_gas.vmr.size();
      if (is_active_(this_gas.molecule, active_gas_list)) {
	// Link to state vector
	this_gas.is_active = true;
	this_gas.molar_abs_conc.clear();
	this_gas.molar_abs_conc >>= x(range(ix,ix+nconc*nt_*np_*ng_-1)).reshape(dimensions(nconc,nt_,np_,ng_));
	LOG << "    ACTIVE: preparing to optimize array of size " << this_gas.molar_abs_conc.dimensions() << "\n";
	if (have_min_max) {
	  this_gas.min_molar_abs_conc >>= x_min(range(ix,ix+nconc*nt_*np_*ng_-1)).reshape(dimensions(nconc,nt_,np_,ng_));
	  this_gas.max_molar_abs_conc >>= x_max(range(ix,ix+nconc*nt_*np_*ng_-1)).reshape(dimensions(nconc,nt_,np_,ng_));
	}
	this_gas.ix = ix;
	ix += nconc*nt_*np_*ng_;
      }
      else {
	this_gas.molar_abs_conc.resize(nconc,nt_,np_,ng_);
	if (have_min_max) {
	  this_gas.min_molar_abs_conc.resize(nconc,nt_,np_,ng_);
	  this_gas.max_molar_abs_conc.resize(nconc,nt_,np_,ng_);
	}
      }
      Array3 tmp;
      for (int iconc = 0; iconc < nconc; ++iconc) {
	file.read(tmp, molecule + "_" + K_NAME, iconc);
	this_gas.molar_abs_conc(iconc,__,__,__) = tmp;
      }
      // If present, read the minimum and maximum possible values of
      // the absorption coefficient
      if (have_min_max) {
	for (int iconc = 0; iconc < nconc; ++iconc) {
	  file.read(tmp, molecule + "_" + K_NAME + "_min", iconc);
	  this_gas.min_molar_abs_conc(iconc,__,__,__) = tmp;
	}
	for (int iconc = 0; iconc < nconc; ++iconc) {
	  file.read(tmp, molecule + "_" + K_NAME + "_max", iconc);
	  this_gas.max_molar_abs_conc(iconc,__,__,__) = tmp;
	}
      }
    }
    else {
      int conc_dep;
      file.read(conc_dep, molecule + "_conc_dependence_code");
      if (conc_dep == 0) {
	this_gas.conc_dependence = NONE;
	LOG << "    Concentration dependence: none\n";
	file.read(this_gas.composite_vmr, molecule + "_mole_fraction");
	file.read(this_gas.composite_molecules, DATA_FILE_GLOBAL_SCOPE,
		  molecule + "_constituent_id");
      }
      else if (conc_dep == 1) {
	this_gas.conc_dependence = LINEAR;
	LOG << "    Concentration dependence: linear\n";
      }
      else if (conc_dep == 3) {
	this_gas.conc_dependence = RELATIVE_LINEAR;
	file.read(this_gas.reference_vmr, molecule + "_reference_mole_fraction");
	LOG << "    Concentration dependence: relative-linear\n";
      }
      else {
	ERROR << molecule + "_conc_dependence_code is inconsistent with other variables in file";
	THROW(1);
      }
      Array3 tmp;
      file.read(tmp, molecule + "_" + K_NAME);

      this_gas.molar_abs.clear();
      if (is_active_(this_gas.molecule, active_gas_list)) {
	this_gas.is_active = true;
	this_gas.molar_abs >>= x(range(ix,ix+nt_*np_*ng_-1)).reshape(dimensions(nt_,np_,ng_));
	LOG << "    ACTIVE: preparing to optimize array of size " << this_gas.molar_abs.dimensions() << "\n";
	if (have_min_max) {
	  this_gas.min_molar_abs >>= x_min(range(ix,ix+nt_*np_*ng_-1)).reshape(dimensions(nt_,np_,ng_));
	  this_gas.max_molar_abs >>= x_max(range(ix,ix+nt_*np_*ng_-1)).reshape(dimensions(nt_,np_,ng_));
	}
	this_gas.ix = ix;
	ix += nt_*np_*ng_;
      }
      this_gas.molar_abs = tmp;

      // If present, read the minimum and maximum possible values of
      // the absorption coefficient
      if (have_min_max) {
	file.read(tmp, molecule + "_" + K_NAME + "_min");
	this_gas.min_molar_abs = tmp;
	file.read(tmp, molecule + "_" + K_NAME + "_max");
	this_gas.max_molar_abs = tmp;
      }
      else {
	this_gas.min_molar_abs.clear();
	this_gas.max_molar_abs.clear();
      }

    }
    ++igas;
  }

  rayleigh_is_active_ = false;
  if (is_sw()) {
    if (file.exist("rayleigh_molar_scattering_coeff")) {
      if (is_active_("rayleigh", active_gas_list)) {
	rayleigh_molar_scat_ >>= x(range(ix,ix+ng_-1));
	rayleigh_is_active_ = true;
	rayleigh_ix_ = ix;
	ix += ng_;
      }
      else {
	rayleigh_molar_scat_.clear();
      }
      Vector tmp;
      file.read(tmp, "rayleigh_molar_scattering_coeff");
      rayleigh_molar_scat_ = tmp;
    }
    else {
      ERROR << "rayleigh_molar_scattering_coeff not present";
      THROW(PARAMETER_ERROR);
    }
  }

  if (!active_gas_list.empty()) {
    active_molecules = active_gas_list;
  }
  else {
    active_molecules = molecules;
  }

  /*
  LOG << "Retrieving gases:";
  for (int igas = 0; igas < active_molecules.size(); ++igas) {
    LOG << " " << active_molecules[igas];
  }
  LOG << "\n";
  */

  if (nx != ix) {
    ERROR << "Mismatch between number of state variables and number of coefficients to optimize";
    THROW(UNEXPECTED_EXCEPTION);
  }

}

template<bool IsActive>
void
CkdModel<IsActive>::write(const std::string& file_name,
			  int argc, const char* argv[], const std::string& config_str)
{
  OutputDataFile file(file_name);

  file.define_dimension("temperature", nt_);
  file.define_dimension("pressure", np_);
  file.define_dimension("g_point", ng_);
  if (!is_sw()) {
    file.define_dimension("temperature_planck", temperature_planck_.size());
  }
  file.define_dimension("wavenumber", nwav_);
  file.define_dimension("band", wavenumber1_band_.size());

  if (do_save_g_points_) {
    file.define_dimension("wavenumber_hr", wavenumber_hr_.size());
  }

  std::string molecule_list;
  for (int igas = 0; igas < ngas(); ++igas) {
    if (igas == 0) {
      molecule_list = single_gas_data_[igas].molecule;
    }
    else {
      molecule_list += " " + single_gas_data_[igas].molecule;
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

  if (is_sw()) {
    file.define_variable("solar_irradiance", FLOAT, "g_point");
    file.write_long_name("Solar irradiance across each g point", "solar_irradiance");
    file.write_units("W m-2", "solar_irradiance");

    if (!ssi_.empty()) {
      file.define_variable("solar_spectral_irradiance", FLOAT, "wavenumber");
      file.write_long_name("Solar irradiance in each spectral interval", "solar_spectral_irradiance");
      file.write_units("W m-2", "solar_spectral_irradiance");
    }
  }
  else {
    file.define_variable("temperature_planck", FLOAT, "temperature_planck");
    file.write_long_name("Temperature for Planck function look-up table", "temperature_planck");
    file.write_units("K", "temperature_planck");

    file.define_variable("planck_function", FLOAT, "temperature_planck", "g_point");
    file.write_long_name("Planck function look-up table", "planck_function");
    file.write_units("W m-2", "planck_function");
  }

  file.define_variable("wavenumber1", FLOAT, "wavenumber");
  file.write_long_name("Lower wavenumber bound of spectral interval", "wavenumber1");
  file.write_units("cm-1", "wavenumber1");
  file.define_variable("wavenumber2", FLOAT, "wavenumber");
  file.write_long_name("Upper wavenumber bound of spectral interval", "wavenumber2");
  file.write_units("cm-1", "wavenumber2");
  file.define_variable("gpoint_fraction", FLOAT, "g_point", "wavenumber");
  file.write_long_name("Fraction of spectrum contributing to each g-point",
		       "gpoint_fraction");

  file.define_variable("wavenumber1_band", FLOAT, "band");
  file.write_long_name("Lower wavenumber bound of band", "wavenumber1_band");
  file.write_units("cm-1", "wavenumber1_band");
  file.define_variable("wavenumber2_band", FLOAT, "band");
  file.write_long_name("Upper wavenumber bound of band", "wavenumber2_band");
  file.write_units("cm-1", "wavenumber2_band");
  file.define_variable("band_number", SHORT, "g_point");
  file.write_long_name("Band number of each g point", "band_number");

  if (do_save_g_points_) {
    file.define_variable("wavenumber_hr", DOUBLE, "wavenumber_hr");
    file.write_long_name("High-resolution wavenumber", "wavenumber_hr");
    file.write_units("cm-1", "wavenumber_hr");
    file.define_variable("g_point", SHORT, "wavenumber_hr");
    file.write_long_name("G point", "g_point");
  }

  if (is_sw()) {
    write_standard_attributes(file, "Definition of a correlated k-distribution model for shortwave gas absorption");

    file.define_variable("rayleigh_molar_scattering_coeff", FLOAT, "g_point");
    file.write_long_name("Rayleigh molar scattering coefficient in each g-point",
			 "rayleigh_molar_scattering_coeff");
    file.write_units("m2 mol-1", "rayleigh_molar_scattering_coeff");
  }
  else {
    write_standard_attributes(file, "Definition of a correlated k-distribution model for longwave gas absorption");
  }
  if (!model_id_.empty()) {
    file.write(model_id_, "model_id");
  }
  file.write(molecule_list, "constituent_id");

  for (int igas = 0; igas < ngas(); ++igas) {
    const SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
    const std::string& Molecule = this_gas.Molecule;
    const std::string& molecule = this_gas.molecule;
    std::string varname = molecule + "_" + K_NAME;

    file.define_variable(molecule + "_conc_dependence_code", SHORT);
    file.write_long_name(Molecule + " concentration dependence code",
			 molecule + "_conc_dependence_code");
    file.write("0: No dependence of absorption on concentration (background gases)\n"
	       "1: Absorption varies linearly with concentration\n"
	       "2: Look-up table for concentration-dependence of absorption\n"
	       "3: Linear dependence on concentration minus a reference value",
	       molecule + "_conc_dependence_code", "definition");

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	file.define_variable(varname, FLOAT, "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of background gases", varname);
	file.write_units("m2 mol-1", varname);
	file.write_comment("This is the absorption cross section of background gases per mole of dry air.",
			   varname);
	if (do_save_min_max_ && !this_gas.min_molar_abs.empty()) {
	  file.define_variable(varname + "_min", FLOAT, "temperature", "pressure", "g_point");
	  file.write_long_name("Minimum molar absorption coefficient of background gases", varname + "_min");
	  file.write_units("m2 mol-1", varname + "_min");
	  file.define_variable(varname + "_max", FLOAT, "temperature", "pressure", "g_point");
	  file.write_long_name("Maximum molar absorption coefficient of background gases", varname + "_max");
	  file.write_units("m2 mol-1", varname + "_max");
	}
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
     case RELATIVE_LINEAR:
      {
	file.define_variable(molecule + "_reference_mole_fraction", FLOAT);
	file.write_long_name("Reference mole fraction of " + Molecule,
			     molecule + "_reference_mole_fraction");
	file.write_units("1", molecule + "_reference_mole_fraction");
	file.write_comment("Subtract this from input mole fractions before multiplying by "
			   + molecule + "_" + K_NAME,
			   molecule + "_reference_mole_fraction");
      }
      // Deliberately fall through to the next item
   case LINEAR:
      {
	file.define_variable(molecule + "_" + K_NAME,
			     FLOAT, "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of " + Molecule,
			     molecule + "_" + K_NAME);
	file.write_units("m2 mol-1", molecule + "_" + K_NAME);
	if (do_save_min_max_ && !this_gas.min_molar_abs.empty()) {
	  file.define_variable(varname + "_min", FLOAT,
			       "temperature", "pressure", "g_point");
	  file.write_long_name("Minimum molar absorption coefficient of " + Molecule,
			       varname + "_min");
	  file.write_units("m2 mol-1", varname + "_min");
	  file.define_variable(varname + "_max", FLOAT,
			       "temperature", "pressure", "g_point");
	  file.write_long_name("Maximum molar absorption coefficient of " + Molecule,
			       varname + "_max");
	  file.write_units("m2 mol-1", varname + "_max");
	}
      }
      break;
    case LUT:
      {
	file.define_dimension(molecule + "_mole_fraction", this_gas.vmr.size());

	file.define_variable(molecule + "_mole_fraction", FLOAT, molecule + "_mole_fraction");
	file.write_long_name(Molecule + " mole fraction for look-up table", molecule + "_mole_fraction");
	file.write_units("1", molecule + "_mole_fraction");

	file.define_variable(molecule + "_" + K_NAME,
			     FLOAT, molecule + "_mole_fraction", "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of " + Molecule,
			     molecule + "_" + K_NAME);
	file.write_units("m2 mol-1", molecule + "_" + K_NAME);
	if (do_save_min_max_ && !this_gas.min_molar_abs_conc.empty()) {
	  file.define_variable(molecule + "_" + K_NAME + "_min", FLOAT,
			       molecule + "_mole_fraction", "temperature", "pressure", "g_point");
	  file.write_long_name("Minimum molar absorption coefficient of " + Molecule,
			       molecule + "_" + K_NAME + "_min");
	  file.write_units("m2 mol-1", molecule + "_" + K_NAME + "_min");
	  file.define_variable(molecule + "_" + K_NAME + "_max", FLOAT,
			       molecule + "_mole_fraction", "temperature", "pressure", "g_point");
	  file.write_long_name("Maximum molar absorption coefficient of " + Molecule,
			       molecule + "_" + K_NAME + "_max");
	  file.write_units("m2 mol-1", molecule + "_" + K_NAME + "_max");
	}
      }
      break;
    }
  }

  // Copy history attribute from earlier file, if present
  if (!history_.empty()) {
    file.write(history_, "history");
  }
  file.append_history(argc, argv);

  if (!config_.empty()) {
    file.write(config_ + "\n" + config_str, "config");
  }
  else {
    file.write(config_str, "config");
  }
  
  /*
  if (!summary_.empty()) {
    std::string summary = summary_ + "\n - Optimized the absorption coefficients of:";
    for (int igas = 0; igas < active_molecules.size(); ++igas) {
      std::string Molecule = active_molecules[igas];
      std::transform(Molecule.begin(), Molecule.end(), Molecule.begin(), ::toupper);
      summary += " " + Molecule;
    }
    file.write(summary, "summary");
  }
  else {
    std::string summary = "This file was created in steps shown by lines of the history and config attributes, corresponding to the following:\n"
      " - Created file by averaging line-by-line absorption coefficients";
    file.write(summary, "summary");
  }
  */

  if (summary_.empty()) {
    std::string xwave;
    if (is_sw()) {
      xwave = "shortwave";
    }
    else {
      xwave = "longwave";
    }
    // !--------!--------!--------!--------!--------!--------!--------!--------!--------!--------
    summary_ = "This file contains the description of a correlated k-distribution model for computing\n"
      + xwave + " gas absorption in the terrestrial atmosphere.  The molar absorption coefficient\n"
      "of each gas and each g point (k term or spectral interval) is implemented as a look-up\n"
      "table versus temperature, pressure, and optionally mole fraction.  The optical depths of\n"
      "each gas should be summed.  The model was created in a multi-step process as described by\n"
      "each line of the history and config global attributes.";
  }
  file.write(summary_, "summary");

  // Write data

  file.end_define_mode();

  file.write(ngas(), "n_gases");
  file.write(eval(exp(log_pressure_)), "pressure");
  file.write(temperature_, "temperature");
  if (is_sw()) {
    file.write(solar_irradiance_, "solar_irradiance");
    file.write(rayleigh_molar_scat_.inactive_link(), "rayleigh_molar_scattering_coeff");
    if (!ssi_.empty()) {
      file.write(ssi_, "solar_spectral_irradiance");
    }
  }   
  else {
    file.write(temperature_planck_, "temperature_planck");
    file.write(planck_function_, "planck_function");
  }
  file.write(wavenumber1_, "wavenumber1");
  file.write(wavenumber2_, "wavenumber2");
  file.write(gpoint_fraction_, "gpoint_fraction");
  file.write(wavenumber1_band_, "wavenumber1_band");
  file.write(wavenumber2_band_, "wavenumber2_band");
  file.write(band_number_, "band_number");

  if (do_save_g_points_) {
    file.write(wavenumber_hr_, "wavenumber_hr");
    file.write(g_point_, "g_point");
  }

  for (int igas = 0; igas < ngas(); ++igas) {
    SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
    const std::string& molecule = this_gas.molecule;

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	file.write(0, molecule + "_conc_dependence_code");
	file.write(this_gas.composite_vmr, molecule + "_mole_fraction");
	file.write((this_gas.molar_abs).inactive_link(), molecule + "_" + K_NAME);
	if (do_save_min_max_ && !this_gas.min_molar_abs.empty()) {
	  file.write((this_gas.min_molar_abs).inactive_link(), molecule + "_" + K_NAME + "_min");
	  file.write((this_gas.max_molar_abs).inactive_link(), molecule + "_" + K_NAME + "_max");
	}
      }
      break;
    case RELATIVE_LINEAR:
      {
	file.write(3, molecule + "_conc_dependence_code");
	file.write(this_gas.reference_vmr, molecule + "_reference_mole_fraction");
	file.write((this_gas.molar_abs).inactive_link(), molecule + "_" + K_NAME);
	if (do_save_min_max_ && !this_gas.min_molar_abs.empty()) {
	  file.write((this_gas.min_molar_abs).inactive_link(), molecule + "_" + K_NAME + "_min");
	  file.write((this_gas.max_molar_abs).inactive_link(), molecule + "_" + K_NAME + "_max");
	}
      }
      break;
    case LINEAR:
      {
	file.write(1, molecule + "_conc_dependence_code");
	file.write((this_gas.molar_abs).inactive_link(), molecule + "_" + K_NAME);
	if (do_save_min_max_ && !this_gas.min_molar_abs.empty()) {
	  file.write((this_gas.min_molar_abs).inactive_link(), molecule + "_" + K_NAME + "_min");
	  file.write((this_gas.max_molar_abs).inactive_link(), molecule + "_" + K_NAME + "_max");
	}
      }
      break;
    case LUT:
      {
	file.write(2, molecule + "_conc_dependence_code");
	file.write(this_gas.vmr, molecule + "_mole_fraction");
	for (int iconc = 0; iconc < this_gas.vmr.size(); ++iconc) {
	  file.write((this_gas.molar_abs_conc(iconc,__,__,__)).inactive_link(),
		     molecule + "_" + K_NAME, iconc);
	  if (do_save_min_max_ && !this_gas.min_molar_abs_conc.empty()) {
	    file.write((this_gas.min_molar_abs_conc(iconc,__,__,__)).inactive_link(),
		       molecule + "_" + K_NAME + "_min", iconc);
	    file.write((this_gas.max_molar_abs_conc(iconc,__,__,__)).inactive_link(),
		       molecule + "_" + K_NAME + "_max", iconc);
	  }
	}
      }
    }
  }

  file.close();

}

/// Create error covariance matrices
template<>
void
CkdModel<true>::create_error_covariances(Real prior_error, Real min_prior_error, 
					 Real max_prior_error, Real prior_error_scaling,
					 Real pressure_corr,
					 Real temperature_corr, Real conc_corr,
					 Real rayleigh_prior_error)
{

  // Exploit sparsity and treat numbers less than this as zero
  static const Real MIN_ERROR_COVARIANCE = 1.0e-6;

  if (prior_error > 0.0) {
    LOG << "Creating a-priori error covariance matrices with fixed error of " << prior_error << "\n";
  }
  else {
    LOG << "Creating a-priori error covariance matrices, estimating errors from minimum and maximum molar absorptions for each g-point\n";
    if (min_prior_error > 0.0) {
      LOG << "  Lower bound on prior error: " << min_prior_error << "\n";
    }
    if (max_prior_error > 0.0) {
      LOG << "  Upper bound on prior error: " << max_prior_error << "\n";
    }
    if (prior_error_scaling > 0.0) {
      LOG << "  Scaling prior error by " << prior_error_scaling << "\n";
    }
  }

  for (int igas = 0; igas < ngas(); ++igas) {
    SingleGasData<true>& this_gas = single_gas_data_[igas];
    if (!this_gas.is_active) {
      continue;
    }

    this_gas.background_error.resize(ng_);
    if (prior_error > 0.0) {
      this_gas.background_error = prior_error;
    }
    else {
      // Default error
      this_gas.background_error = 1.0;
    }

    if (this_gas.conc_dependence == LINEAR
	|| this_gas.conc_dependence == NONE
	|| this_gas.conc_dependence == RELATIVE_LINEAR) {
      int nx = nt_ * np_;
      LOG << "  Creating " << nx << "x" << nx << " error covariance matrix for " << this_gas.Molecule << "\n";

      intVector t_index_vec(nx), p_index_vec(nx);
      // Soft links:
      intMatrix t_index = t_index_vec.reshape(dimensions(nt_,np_));
      intMatrix p_index = p_index_vec.reshape(dimensions(nt_,np_));

      // SPREAD CONVENTION IS WRONG!!!
      //      t_index = spread<0>(range(0,nt_-1),np_);
      //      p_index = spread<1>(range(0,np_-1),nt_);
      t_index = spread<1>(range(0,nt_-1),np_);
      p_index = spread<0>(range(0,np_-1),nt_);

      SymmMatrix background(nx,nx);

      background = pow(temperature_corr,1.0*abs(spread<0>(t_index_vec,nx)-spread<1>(t_index_vec,nx)));
      background *= pow(pressure_corr,1.0*abs(spread<0>(p_index_vec,nx)-spread<1>(p_index_vec,nx)));

      Matrix inv_background = inv(background);
      LOG << "    fraction of elements less than " << MIN_ERROR_COVARIANCE << " is "
	  << static_cast<Real>(count(fabs(inv_background) < MIN_ERROR_COVARIANCE))
	/static_cast<Real>(inv_background.size()) << "\n";
      inv_background.where(fabs(inv_background) < MIN_ERROR_COVARIANCE) = 0.0;
      // Save as a sparse matrix (if supported in later versions of
      // Adept) or copy to a symmetric matrix
      this_gas.inv_background_shape = inv_background;

      // Estimate the prior error of the logarithm of absorption for
      // each g point as 0.25*[log(max_molar_abs)-log(min_molar_abs)],
      // or if min is zero, 0.5*[log(max_molar_abs)-log(molar_abs)]
      if (prior_error <= 0.0) {
	for (int ig = 0; ig < ng_; ++ig) {
	  Real local_sum = 0.0;
	  int local_count = 0;
	  for (int ip = 0; ip < np_; ++ip) {
	    for (int it = 0; it < nt_; ++it) {
	      if (this_gas.molar_abs(it,ip,ig) > 0.0) {
		if (this_gas.min_molar_abs(it,ip,ig) > 0.0) {
		  local_sum += 0.25 * log(this_gas.max_molar_abs(it,ip,ig)/this_gas.min_molar_abs(it,ip,ig));
		}
		else {
		  local_sum += 0.5 * log(this_gas.max_molar_abs(it,ip,ig)/value(this_gas.molar_abs(it,ip,ig)));
		}
		++local_count;
	      }
	    }
	  }
	  if (local_count > 0) {
	    this_gas.background_error(ig) = prior_error_scaling * local_sum / local_count;
	  }
	}
	if (min_prior_error > 0.0) {
	  this_gas.background_error = max(min_prior_error, this_gas.background_error);
	}
	if (max_prior_error > 0.0) {
	  this_gas.background_error = min(this_gas.background_error, max_prior_error);
	}
	LOG << "    prior errors for " << this_gas.Molecule << " per g-point: "
	    << this_gas.background_error<< "\n";
      }
    }
    else {
      int nconc = this_gas.vmr.size();
      int nx = nt_ * np_ * nconc;
      LOG << "  Creating " << nx << "x" << nx << " error covariance matrix for " << this_gas.Molecule << "\n";
      intVector t_index_vec(nx), p_index_vec(nx), c_index_vec(nx);
      // Soft links:
      intArray3D t_index = t_index_vec.reshape(dimensions(nconc,nt_,np_));
      intArray3D p_index = p_index_vec.reshape(dimensions(nconc,nt_,np_));
      intArray3D c_index = c_index_vec.reshape(dimensions(nconc,nt_,np_));

      // SPREAD CONVENTION IS WRONG!!!
      //      t_index = spread<2>(spread<0>(range(0,nt_-1),np_),nconc);
      //      p_index = spread<2>(spread<1>(range(0,np_-1),nt_),nconc);
      //      c_index = spread<0>(spread<0>(range(0,nconc-1),nt_),np_);
      t_index = spread<0>(spread<1>(range(0,nt_-1),np_),nconc);
      p_index = spread<0>(spread<0>(range(0,np_-1),nt_),nconc);
      c_index = spread<2>(spread<1>(range(0,nconc-1),nt_),np_);

      Matrix background(nx,nx);

      background = pow(temperature_corr,1.0*abs(spread<0>(t_index_vec,nx)-spread<1>(t_index_vec,nx)));
      background *= pow(pressure_corr,1.0*abs(spread<0>(p_index_vec,nx)-spread<1>(p_index_vec,nx)));
      background *= pow(conc_corr,1.0*abs(spread<0>(c_index_vec,nx)-spread<1>(c_index_vec,nx)));

      Matrix inv_background = inv(background);
      LOG << "    fraction of elements less than " << MIN_ERROR_COVARIANCE << " is "
	  << static_cast<Real>(count(fabs(inv_background) < MIN_ERROR_COVARIANCE))
	/static_cast<Real>(inv_background.size()) << "\n";
      inv_background.where(fabs(inv_background) < MIN_ERROR_COVARIANCE) = 0.0;
      // Save as a sparse matrix
      this_gas.inv_background_shape = inv_background;

      // Estimate the prior error of the logarithm of absorption for
      // each g point as 0.25*[log(max_molar_abs)-log(min_molar_abs)],
      // or if min is zero, 0.5*[log(max_molar_abs)-log(molar_abs)]
      if (prior_error <= 0.0) {
	int nconc = this_gas.vmr.size();
	for (int ig = 0; ig < ng_; ++ig) {
	  Real local_sum = 0.0;
	  int local_count = 0;
	  for (int ic = 0; ic < nconc; ++ic) {
	    for (int ip = 0; ip < np_; ++ip) {
	      for (int it = 0; it < nt_; ++it) {
		if (this_gas.molar_abs_conc(ic,it,ip,ig) > 0.0) {
		  if (this_gas.min_molar_abs_conc(ic,it,ip,ig) > 0.0) {
		    local_sum += 0.25 * log(this_gas.max_molar_abs_conc(ic,it,ip,ig)
					   /this_gas.min_molar_abs_conc(ic,it,ip,ig));
		  }
		  else {
		    local_sum += 0.5 * log(this_gas.max_molar_abs_conc(ic,it,ip,ig)
					   /value(this_gas.molar_abs_conc(ic,it,ip,ig)));
		  }
		  ++local_count;
		}
	      }
	    }
	  }
	  if (local_count > 0) {
	    this_gas.background_error(ig) = prior_error_scaling * local_sum / local_count;
	  }
	}
	if (min_prior_error > 0.0) {
	  this_gas.background_error = max(min_prior_error, this_gas.background_error);
	}
	if (max_prior_error > 0.0) {
	  this_gas.background_error = min(this_gas.background_error, max_prior_error);
	}
	LOG << "    prior errors for " << this_gas.Molecule << " per g-point: "
	    << this_gas.background_error<< "\n";
     }
    }
  }
  if (rayleigh_prior_error > 0.0 && rayleigh_is_active_) {
    rayleigh_inv_background_.resize(ng_);
    rayleigh_inv_background_ = 1.0 / (rayleigh_prior_error*rayleigh_prior_error);
  }
  else {
    rayleigh_inv_background_.clear();
  }

}


/// Return the background contribution to cost function, J, and also
/// the gradient dJ/dx, where delta_x is the difference between the
/// current state and the prior
template<>
Real
CkdModel<true>::calc_background_cost_function(const Vector& delta_x, Vector gradient)
{
  Real cost_fn = 0.0;
  gradient = 0.0;

  for (int igas = 0; igas < ngas(); ++igas) {
    const SingleGasData<true>& this_gas = single_gas_data_[igas];
    if (this_gas.is_active) {
      int nx = this_gas.inv_background_shape.dimension(0);
      // g-point is fastest varying dimension so need to stride over it
      int nstride = ng_;
      //#pragma omp parallel for schedule (static)
      for (int ig = 0; ig < ng_; ++ig) {
	int ix = this_gas.ix + ig;
	Vector delta_x_local = delta_x.soft_link()(stride(ix,ix+(nx-1)*nstride,nstride));
	//	LOG << "??? " << this_gas.Molecule << " " << delta_x.size() << " " << ix << " " << nx << " " << nstride << " " << ig << " " << delta_x_local.size() << "\n";
	Vector gradient_local = (1.0/(this_gas.background_error(ig)*this_gas.background_error(ig)))
	  * (this_gas.inv_background_shape.soft_link() ** delta_x_local);
	Real cost_fn_local = 0.5*dot_product(delta_x_local,gradient_local);

	gradient.soft_link()(stride(ix,ix+(nx-1)*nstride,nstride)) += gradient_local;
	//#pragma omp critical
	{
	  cost_fn += cost_fn_local;
	}
      }
    }
  }

  if (rayleigh_is_active_ && !rayleigh_inv_background_.empty()) {
    Vector gradient_local(ng_);
    Vector delta_x_local = delta_x(range(rayleigh_ix_,rayleigh_ix_+ng_-1));
    gradient_local = rayleigh_inv_background_ * delta_x_local;
    cost_fn += 0.5*dot_product(delta_x_local,gradient_local);
  }

  return cost_fn;
}

/// Ensure that gases with a "relative-linear" representation cannot
/// lead to a negative optical depth if their concentration is zero
template<>
void
CkdModel<true>::cap_relative_linear_coeffts(Real ref_frac_trigger)
{
  // First find index of background gas
  int ibggas = -1;
  bool is_rel_lin = false;
  for (int igas = 0; igas < ngas(); ++igas) {
    SingleGasData<true>& this_gas = single_gas_data_[igas];
    if (this_gas.conc_dependence == NONE) {
      ibggas = igas;
    }
    else if (this_gas.is_active
	     && this_gas.conc_dependence == RELATIVE_LINEAR) {
      is_rel_lin = true;
    }
  }
  if (is_rel_lin) {
    if (ibggas < 0) {
      LOG << "Unable to cap the coefficients of the relative-linear gases because no background composite gas found\n";
    }
    else {
      for (int igas = 0; igas < ngas(); ++igas) {
	SingleGasData<true>& this_gas = single_gas_data_[igas];
	if (this_gas.is_active
	    && this_gas.conc_dependence == RELATIVE_LINEAR) {
	  int nbad = count(this_gas.molar_abs*(this_gas.reference_vmr*ref_frac_trigger) > single_gas_data_[ibggas].molar_abs);
	  if (nbad > 0) {
	    LOG << "Correcting " << nbad << " " << this_gas.Molecule << " coefficients that could cause negative optical depth\n";
	    this_gas.molar_abs = min(this_gas.molar_abs, single_gas_data_[ibggas].molar_abs/(this_gas.reference_vmr*ref_frac_trigger));
	  }
	}
      }  
    }
  }

}


/// Calculate the optical depth at each g-point for multiple
/// atmospheric columns, where the arguments are dimensioned
/// (column,level) and the output dimensioned (column,level,g-point)
template<bool IsActive>
Array<3,Real,IsActive>
CkdModel<IsActive>::calc_optical_depth(int igas,                         ///< Gas number 
				       const Matrix& pressure_hl,        ///< Pressure at half levels (Pa)
				       const Matrix& temperature_fl,     ///< Temperature at full levels (K)
				       const Matrix& vmr_fl) ///< Volume mixing ratio at full levels
{ 
  using std::fmin;
  using std::fmax;
  
  typedef typename scalar<IsActive>::type areal;
  typedef Array<1,Real,IsActive> avector;
  typedef Array<2,Real,IsActive> amatrix;
  typedef Array<3,Real,IsActive> aarray3D;

  /*
  LOG << molecules[igas] << " " << pressure_hl(0,20) << " " << temperature_fl(0,20);
  if (!vmr_fl.empty()) {
    LOG << " " << vmr_fl(0,20) << "\n";
  }
  else {
    LOG << "\n";
  }
  */

  // Assume pressure of LUT is evenly spaced in log space
  Real log_p_0 = log_pressure_(0);
  Real d_log_p = log_pressure_(1)-log_pressure_(0); // Spacing

  // Temperature LUT spacing
  Real d_t = temperature_(1,0)-temperature_(0,0);

  SingleGasData<IsActive>& this_gas = single_gas_data_[igas];

  int ncol = pressure_hl.dimension(0);
  int np   = pressure_hl.dimension(1)-1;
  aarray3D od(ncol,np,ng_);
  Real global_weight = 1.0 / (ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR);

  for (int icol = 0; icol < ncol; ++icol) {
    for (int ip = 0; ip < np; ++ip) {
      // Find interpolation points in pressure
      Real log_pressure_fl = log(0.5*(pressure_hl(icol,ip+1)+pressure_hl(icol,ip)));
      Real pindex0 = (log_pressure_fl-log_p_0) / d_log_p;
      pindex0 = fmax(0.0, fmin(pindex0, np_-1.0001));
      int ip0 = static_cast<int>(pindex0);
      Real pweight1 = pindex0 - ip0;
      Real pweight0 = 1.0 - pweight1;
      // Find interpolation points in temperature
      Real t_0 = pweight0*temperature_(0,ip0) + pweight1*temperature_(0,ip0+1);
      Real tindex0 = (temperature_fl(icol,ip)-t_0) / d_t;
      tindex0 = fmax(0.0, fmin(tindex0, nt_-1.0001));
      int it0 = static_cast<int>(tindex0);
      Real tweight1 = tindex0 - it0;
      Real tweight0 = 1.0 - tweight1;

      // Weight
      Real simple_weight = global_weight * (pressure_hl(icol,ip+1)-pressure_hl(icol,ip));

      /*
      if (icol == 0) {
	std::cout << ip << " " << log_pressure_fl << " " << pindex0 << " " << ip0 << " " << pweight0 << " " << pweight1 << " " << t_0 << " " << tindex0 << " " << it0 << " " << tweight0 << " " << tweight1 << " " << simple_weight << "\n";
	std::cout << "  " << temperature_fl(icol,ip) << " " << d_t << " " << nt_ << "\n";
      }
      */

      Real weight = 0.0;
      bool no_vmr_provided = true;
      if (!vmr_fl.empty()) {
	if (this_gas.conc_dependence == RELATIVE_LINEAR) {
	  weight = simple_weight * (vmr_fl(icol,ip) - this_gas.reference_vmr);
	}
	else {
	  weight = simple_weight * vmr_fl(icol,ip);
	}
	no_vmr_provided = false;
      }

      if (this_gas.conc_dependence == LUT) {
	// Find interpolation points in concentration
	Real log_conc = log(vmr_fl(icol,ip));
	Real d_log_c  = log(this_gas.vmr(1)/this_gas.vmr(0));
	Real cindex0  = (log_conc-log(this_gas.vmr(0))) / d_log_c;
	cindex0 = fmax(0.0, fmin(cindex0, this_gas.vmr.size()-1.0001));
	int ic0 = static_cast<int>(cindex0);
	Real cweight1 = cindex0 - ic0;
	Real cweight0 = 1.0 - cweight1;
	/*
        if (icol == 0) {
	  std::cout << "  " << log_conc << " " << vmr_fl(icol,ip) << " " << cindex0 << " " << log(this_gas.vmr(0)) << " " << ic0 << " " << cweight0 << " " << cweight1 << " " << d_log_c << "\n";
	}
	*/

	if (no_vmr_provided) {
	  ERROR << "Concentration of " << molecules[igas] << " not provided";
	  THROW(1);
	}
	// Tri-linear interpolation
	if (!logarithmic_interpolation) {
	  od(icol,ip,__) = weight
	    * (cweight0
	       * (tweight0  * (pweight0 * this_gas.molar_abs_conc(ic0,it0,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0,it0,ip0+1,__))
		  +tweight1 * (pweight0 * this_gas.molar_abs_conc(ic0,it0+1,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0,it0+1,ip0+1,__)))
	       +cweight1
	       * (tweight0  * (pweight0 * this_gas.molar_abs_conc(ic0+1,it0,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0+1,it0,ip0+1,__))
		  +tweight1 * (pweight0 * this_gas.molar_abs_conc(ic0+1,it0+1,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0+1,it0+1,ip0+1,__))));
	}
	else {
	  od(icol,ip,__) = weight
	    * exp(cweight0
		  * (tweight0  * (pweight0 * log(this_gas.molar_abs_conc(ic0,it0,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0,it0,ip0+1,__)))
		     +tweight1 * (pweight0 * log(this_gas.molar_abs_conc(ic0,it0+1,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0,it0+1,ip0+1,__))))
		  +cweight1
		  * (tweight0  * (pweight0 * log(this_gas.molar_abs_conc(ic0+1,it0,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0+1,it0,ip0+1,__)))
		     +tweight1 * (pweight0 * log(this_gas.molar_abs_conc(ic0+1,it0+1,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0+1,it0+1,ip0+1,__)))));
	}
      }
      else if (this_gas.conc_dependence == LINEAR
	       || this_gas.conc_dependence == RELATIVE_LINEAR) {
	if (no_vmr_provided) {
	  ERROR << "Concentration of " << molecules[igas] << " not provided";
	  THROW(1);
	}
	if (!logarithmic_interpolation) {
	  // Bi-linear interpolation
	  od(icol,ip,__) = weight
	    * (tweight0  * (pweight0 * this_gas.molar_abs(it0,ip0,__)
			    +pweight1* this_gas.molar_abs(it0,ip0+1,__))
	       +tweight1 * (pweight0 * this_gas.molar_abs(it0+1,ip0,__)
			    +pweight1* this_gas.molar_abs(it0+1,ip0+1,__)));
	}
	else {
	  od(icol,ip,__) = weight
	    * exp(tweight0  * (pweight0 * log(this_gas.molar_abs(it0,ip0,__))
			       +pweight1* log(this_gas.molar_abs(it0,ip0+1,__)))
		  +tweight1 * (pweight0 * log(this_gas.molar_abs(it0+1,ip0,__))
			       +pweight1* log(this_gas.molar_abs(it0+1,ip0+1,__))));
	}
      }
      else { // NONE
	if (!logarithmic_interpolation) {
	  // Bi-linear interpolation
	  od(icol,ip,__) = simple_weight
	    * (tweight0  * (pweight0 * this_gas.molar_abs(it0,ip0,__)
			    +pweight1* this_gas.molar_abs(it0,ip0+1,__))
	       +tweight1 * (pweight0 * this_gas.molar_abs(it0+1,ip0,__)
			    +pweight1* this_gas.molar_abs(it0+1,ip0+1,__)));
	}
	else {
	  od(icol,ip,__) = simple_weight
	    * exp(tweight0  * (pweight0 * log(this_gas.molar_abs(it0,ip0,__))
			       +pweight1* log(this_gas.molar_abs(it0,ip0+1,__)))
		  +tweight1 * (pweight0 * log(this_gas.molar_abs(it0+1,ip0,__))
			       +pweight1* log(this_gas.molar_abs(it0+1,ip0+1,__))));
	}
      }
    }
  }

  /*
  LOG << "  CKD_MODEL: requested optical depth of " << this_gas.molecule
      << " " << pressure_hl(0,end) << " " << temperature_fl(0,end);
  if (!vmr_fl.empty()) {
    LOG << " [" << vmr_fl(0,end) << "]";
  }

  LOG << " " << od(0,end,0) << "\n";
  */
  //  LOG << "abs " << molecules[igas] << " = " << od(0,20,__) << "\n";

  return od;
}


template<bool IsActive>
Array3
CkdModel<IsActive>::calc_planck_function(const Matrix& temperature_hl)
{
  int ncol = temperature_hl.dimension(0);
  int np   = temperature_hl.dimension(1);
  Array3 planck(ncol,np,ng_);
  for (int icol = 0; icol < ncol; ++icol) {
    planck(icol,__,__) = calc_planck_function(temperature_hl(icol,__));
  }
  return planck;
}


template<bool IsActive>
Matrix
CkdModel<IsActive>::calc_planck_function(const Vector& temperature)
{
  // Planck temperature LUT spacing
  Real d_t = temperature_planck_(1) - temperature_planck_(0);
  Real t0 = temperature_planck_(0);
  int nt = temperature.size();
  Matrix planck(nt,ng_);
  for (int it = 0; it < nt; ++it) {
    Real tindex0 = (temperature(it)-t0) / d_t;
    if (tindex0 >= 0) {
      // Normal interpolation (and extrapolation for high
      // temperatures)
      int it0 = std::min(static_cast<int>(tindex0), temperature_planck_.size()-2);
      Real tweight1 = tindex0 - it0;
      Real tweight0 = 1.0 - tweight1;
      planck(it,__) = tweight0 * planck_function_(it0,__)
   	            + tweight1 * planck_function_(it0+1,__);
    }
    else {
      // Interpolate linearly to zero
      planck(it,__) = (temperature(it)/t0) * planck_function_(0,__);
    }
  }
  return planck;
}

/// Scale the optical depth coefficients of each gas equally, where
/// scaling is dimensioned (nz,ng)
template<bool IsActive>
void
CkdModel<IsActive>::scale_optical_depth(const Vector& pressure_fl, const Matrix& scaling) {
  Matrix local_scaling = interp(log(pressure_fl), scaling, log_pressure_);
  for (int ip = 2; ip < log_pressure_.size(); ip += 10) {
    LOG << "  Scalings at " << exp(log_pressure_(ip)) << " Pa: " << local_scaling[ip] << "\n";
  }

  for (int igas = 0; igas < ngas(); ++igas) {
    SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
    if (this_gas.conc_dependence == LUT) {
      this_gas.molar_abs_conc *= spread<0>(spread<0>(local_scaling, nt_), this_gas.molar_abs_conc.size(0));
      // If minimum and maximum possible absorption coefficients are
      // available, check we are still within these bounds
      if (!this_gas.min_molar_abs_conc.empty()) {
	this_gas.molar_abs_conc = max(this_gas.min_molar_abs_conc, 
				      min(this_gas.molar_abs_conc, this_gas.max_molar_abs_conc));
      }
    }
    else {
      this_gas.molar_abs *= spread<0>(local_scaling, nt_);
      if (!this_gas.min_molar_abs.empty()) {
	this_gas.molar_abs = max(this_gas.min_molar_abs, 
				 min(this_gas.molar_abs, this_gas.max_molar_abs));
      }
    }
  }
}


// Instantiate both active and passive versions of CkdModel

template CkdModel<false>::CkdModel(const std::vector<SingleGasData<false> >&,
				   const Vector&, const Matrix&, const Vector&, const Matrix&,
				   const Vector&, const Vector&, const Matrix&,
				   const Vector& wavenumber1_band,
				   const Vector& wavenumber2_band,
				   const intVector& band_number,
				   const std::string& history,
				   const std::string& config);
//template CkdModel<true>::CkdModel(const std::vector<SingleGasData<false> >&,
//				  const Vector&, const Matrix&, const Vector&, const Matrix&,
//				  const Vector&, const Vector&, const Matrix&);
template void CkdModel<false>::read(const std::string&,const std::vector<std::string>&);
template void CkdModel<true>::read(const std::string&,const std::vector<std::string>&);
template void CkdModel<false>::write(const std::string&, int argc, const char* argv[], const std::string&);
template void CkdModel<true>::write(const std::string&, int argc, const char* argv[], const std::string&);
template Array<3,Real,false> CkdModel<false>::calc_optical_depth(int,const Matrix&,const Matrix&,const Matrix&);
template Array<3,Real,true> CkdModel<true>::calc_optical_depth(int,const Matrix&,const Matrix&,const Matrix&);
template Array3 CkdModel<false>::calc_planck_function(const Matrix&);
template Array3 CkdModel<true>::calc_planck_function(const Matrix&);
template Matrix CkdModel<false>::calc_planck_function(const Vector&);
template Matrix CkdModel<true>::calc_planck_function(const Vector&);
template void CkdModel<false>::scale_optical_depth(const Vector&, const Matrix&);
