// read_merged_spectrum.cpp - Read and combine the spectral optical depths of several gases
//
// Copyright (C) 2019- ECMWF.
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

#include "read_merged_spectrum.h"
#include "read_spectrum.h"

/// Read and combine the spectral optical depths of several gases
void
read_merged_spectrum(DataFile& config,                ///< Config file
		     int iprofile,                    ///< Index of profile (0 based)
		     std::string prefix,              ///< Prefix of keys used to query "config"
		     adept::Vector& pressure_hl,      ///< Half-level pressure (Pa)
		     adept::Vector& temperature_hl,   ///< Half-level temperature (K)
		     adept::Vector& wavenumber_cm_1,  ///< Wavenumber (cm-1)
		     adept::Vector& d_wavenumber_cm_1,///< Wavenumber interval (cm-1)
		     adept::Matrix& optical_depth,    ///< Spectral optical depth 
		     std::string& molecules,          ///< Molecule formulas separated by hyphens
		     adept::Matrix& vmr_fl,           ///< Volume mixing ratio (mol mol-1)
		     int* ngas,                       ///< Number of gases read in
		     int* ncol                        ///< Number of columns in file
		     ) {

  // Names of keys to look for in config file
  std::string input_name   = prefix + "input";
  std::string scaling_name = prefix + "scaling";
  std::string conc_name    = prefix + "conc";
  std::string conc_input_name = prefix + "conc_input";
  std::string iprof_conc_name = prefix + "iprofile";

  std::string file_name;

  int ibg = 0; // Gas index

  intVector tmpvec = config.size(input_name);
  int num_gases = tmpvec(0);

  int iprof_conc = -1;

  // Do we have the name of a concentration file to be used to scale
  // the spectra?
  std::string conc_file_name;
  DataFile conc_file;
  Vector pressure_conc;
  if (config.read(conc_file_name, conc_input_name)) {
    if (!config.read(iprof_conc, iprof_conc_name)) {
      ERROR << "Concentration file specified without profile number in \"iprofile\"";
      THROW(PARAMETER_ERROR);
    }
    conc_file.open(conc_file_name);
    conc_file.read(pressure_conc,"pressure_fl",iprof_conc);
  }

  Vector pressure_fl;
  while(config.read(file_name, input_name, ibg)) {
    Real scaling = -1.0;
    Real conc = -1.0;
    config.read(scaling, scaling_name, ibg);
    config.read(conc, conc_name, ibg);
    Matrix od;
    Real reference_surface_vmr;
    Vector vmr_fl_one_gas;
    std::string molecule;
    Vector scaling_profile;
    Vector conc_interp;

    LOG << "  Reading " << file_name << "\n";
    if (ibg == 0) {
      // First spectrum
      read_spectrum(file_name, iprofile, pressure_hl, temperature_hl,
		    wavenumber_cm_1, d_wavenumber_cm_1, od, 
		    molecule, reference_surface_vmr, vmr_fl_one_gas, ncol);
      molecules = molecule;
      pressure_fl.resize(pressure_hl.size()-1);
      pressure_fl = 0.5*(pressure_hl(range(0,end-1)) + pressure_hl(range(1,end)));
    }
    else {
      // Subsequent spectra
      DataFile file(file_name);
      file.read(od, "optical_depth", iprofile);
      if (file.exist("reference_surface_mole_fraction")) {
	file.read(reference_surface_vmr, "reference_surface_mole_fraction");
      }
      else {
	reference_surface_vmr = -1.0;
      }
      // Read volume mixing ratio on full levels, or negative value if not
      // present (e.g. for hybrid spectra)
      if (file.exist("mole_fraction_fl") && file.size("mole_fraction_fl").size() == 2) {
	file.read(vmr_fl_one_gas, "mole_fraction_fl", iprofile);
      }
      else {
	vmr_fl_one_gas.resize(pressure_hl.size()-1);
	vmr_fl_one_gas = -1.0;
      }
      molecule.clear();
      if (!file.read(molecule, DATA_FILE_GLOBAL_SCOPE, "constituent_id")) {
	file.read(molecule, DATA_FILE_GLOBAL_SCOPE, "molecules");
      }
      if (molecule.empty()) {
	ERROR << "Found neither \"constituent_id\" nor \"molecules\" amongst the global attributes";
	THROW(PARAMETER_ERROR);
      }
      molecules += " " + molecule;
      file.close();
    }
    if (iprof_conc >= 0) {
      // We have a requested concentration profile so need to scale
      // the optical depths of the gases
      Vector conc_req;
      conc_file.read(conc_req, molecule + "_mole_fraction_fl", iprof_conc);
      conc_interp = interp(pressure_conc, conc_req, pressure_fl);
      // Clamp the ends
      conc_interp.where(pressure_fl < pressure_conc(0)) = conc_req(0);
      conc_interp.where(pressure_fl > pressure_conc(end)) = conc_req(end);
      scaling_profile = conc_interp / vmr_fl_one_gas;
      //LOG << "    conc_interp = " << conc_interp << "\n";
      //      LOG << "    scaling_profile = " << scaling_profile << "\n";
      LOG << "    Scaling to target concentration profile in the range "
	  << minval(conc_interp) << " to " << maxval(conc_interp) << "\n";
    }
    else if (conc == 0.0) {
      scaling = 0.0;
    }
    else if (conc > 0.0) {
      if (reference_surface_vmr < 0.0) {
	ERROR << "Attempt to specify concentration when no reference_surface_mole_fraction present in "
	      << file_name;
	THROW(PARAMETER_ERROR);
      }
      scaling = conc / reference_surface_vmr;
      LOG << "    Reference surface concentration = " << reference_surface_vmr << "\n";
      LOG << "    Target surface concentration    = " << conc << "\n";
    }
    else if (scaling < 0.0) {
      scaling = 1.0;
    }

    if (ibg == 0) {
      // First spectrum
      optical_depth.resize(od.dimensions());
      optical_depth = 0.0;
      vmr_fl.resize(num_gases, pressure_hl.size()-1);
    }

    if (iprof_conc >= 0) {
      optical_depth += od * spread<1>(scaling_profile,od.size(1));
      vmr_fl(ibg,__) = conc_interp;
    }
    else if (scaling != 1.0) {
      LOG << "    Scaling by " << scaling << "\n";
      optical_depth += od * scaling;
      vmr_fl(ibg,__) = vmr_fl_one_gas * scaling;
    }
    else {
      optical_depth += od;
      vmr_fl(ibg,__) = vmr_fl_one_gas;
    }
    od.clear();
    ++ibg;
  }

  if (ibg == 0) {
    ERROR << "Unable to read input file names in " << input_name << "\n";
    THROW(PARAMETER_ERROR);
  }

  if (ngas) {
    *ngas = ibg;
  }

  {
    Vector col_od = sum(optical_depth,1);
    Real mean_od = mean(col_od);
    Real mean_od2 = mean(col_od*col_od);
    Real std_od = sqrt(mean_od2 - mean_od*mean_od);
    LOG << "    Column optical depth: " << mean_od << " +/- " << std_od << "\n";
  }

}
