#include <string>
#include <adept_arrays.h>
#include "DataFile.h"

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
	      adept::Real& reference_surface_vmr,///< Reference volume mixing ratio (or -1.0)
	      adept::Vector& vmr_fl,           ///< Volume mixing ratio on full levels
	      int* ncol = 0                    ///< Number of columns in file
	      );
