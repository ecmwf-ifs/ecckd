#include "Error.h"
#include "read_spectrum.h"

/// Read a profile of spectral optical depth from a NetCDF file
void
read_spectrum(std::string& file_name,          ///< File name containing spectra
	      int iprof,                       ///< Index of profile (0 based)
	      adept::Vector& pressure_hl,      ///< Half-level pressure (Pa)
	      adept::Vector& temperature_hl,   ///< Half-level temperature (K)
	      adept::Vector& wavenumber_cm_1,  ///< Wavenumber (cm-1)
	      adept::Vector& d_wavenumber_cm_1,///< Wavenumber interval (cm-1)
	      adept::Matrix& optical_depth,    ///< Spectral optical depth 
	      std::string& molecule,           ///< Chemical formula of molecule 
	      Real& reference_surface_vmr,     ///< Reference volume mixing ratio (or -1.0)
	      adept::Vector& vmr_fl,           ///< Volume mixing ratio on full levels
	      int* ncol                        ///< Number of columns in file
	      ) {

  DataFile file(file_name);

  pressure_hl.clear();
  wavenumber_cm_1.clear();
  d_wavenumber_cm_1.clear();
  optical_depth.clear();

  if (ncol) {
    intVector dims = file.size("pressure_hl");
    *ncol = dims(0);
  }

  file.read(pressure_hl, "pressure_hl", iprof);
  if (file.exist("temperature_hl")) {
    file.read(temperature_hl, "temperature_hl", iprof);
  }
  else {
    WARNING << "\"temperature_hl\" not present";
    ENDWARNING;
  }
  file.read(wavenumber_cm_1, "wavenumber");

  if (file.exist("d_wavenumber")) {
    file.read(d_wavenumber_cm_1, "d_wavenumber");
  }
  else {
    d_wavenumber_cm_1.resize(wavenumber_cm_1.size());
    d_wavenumber_cm_1(range(1,end-1))
      = 0.5 * (wavenumber_cm_1(range(2,end))
	       -wavenumber_cm_1(range(0,end-2)));
    d_wavenumber_cm_1(0) = 0.5*d_wavenumber_cm_1(1);
    d_wavenumber_cm_1(end) = 0.5*d_wavenumber_cm_1(end-1);
  }

  file.read(molecule, DATA_FILE_GLOBAL_SCOPE, "constituent_id");
  if (file.exist("reference_surface_mole_fraction")) {
    file.read(reference_surface_vmr, "reference_surface_mole_fraction");
  }
  else {
    reference_surface_vmr = -1.0;
  }
  // Read volume mixing ratio on full levels, or negative value if not
  // present (e.g. for hybrid spectra)
  if (file.exist("mole_fraction_fl") && file.size("mole_fraction_fl").size() == 2) {
    file.read(vmr_fl, "mole_fraction_fl", iprof);
  }
  else {
    vmr_fl.resize(pressure_hl.size()-1);
    vmr_fl = -1.0;
  }

  file.read(optical_depth, "optical_depth", iprof);

  file.close();
}
