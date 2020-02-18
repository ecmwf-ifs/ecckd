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

  std::string input_name   = prefix + "input";
  std::string scaling_name = prefix + "scaling";
  std::string conc_name    = prefix + "conc";

  std::string file_name;

  int ibg = 0;

  intVector tmpvec = config.size(input_name);
  int num_gases = tmpvec(0);

  while(config.read(file_name, input_name, ibg)) {
    Real scaling = -1.0;
    Real conc = -1.0;
    config.read(scaling, scaling_name, ibg);
    config.read(conc, conc_name, ibg);
    Matrix od;
    Real reference_surface_vmr;
    Vector vmr_fl_one_gas;
    LOG << "  Reading " << file_name << "\n";
    if (ibg == 0) {
      // First spectrum
      read_spectrum(file_name, iprofile, pressure_hl, temperature_hl,
		    wavenumber_cm_1, d_wavenumber_cm_1, od, 
		    molecules, reference_surface_vmr, vmr_fl_one_gas, ncol);
    }
    else {
      // Subsequent spectra
      DataFile file(file_name);
      file.read(od, "optical_depth", iprofile);
      if (file.exist("reference_surface_vmr")) {
	file.read(reference_surface_vmr, "reference_surface_vmr");
      }
      else {
	reference_surface_vmr = -1.0;
      }
      // Read volume mixing ratio on full levels, or negative value if not
      // present (e.g. for hybrid spectra)
      if (file.exist("vmr_fl")) {
	file.read(vmr_fl_one_gas, "vmr_fl", iprofile);
      }
      else {
	vmr_fl_one_gas.resize(pressure_hl.size()-1);
	vmr_fl_one_gas = -1.0;
      }
      std::string mol;
      file.read(mol, DATA_FILE_GLOBAL_SCOPE, "molecule");
      molecules += "," + mol;
      file.close();
    }
    if (conc == 0.0) {
      scaling = 0.0;
    }
    else if (conc > 0.0) {
      if (reference_surface_vmr < 0.0) {
	ERROR << "Attempt to specify concentration when no reference_surface_vmr present in "
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

    if (scaling != 1.0) {
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
}
