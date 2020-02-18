#include <adept_arrays.h>

/// Compute the Planck function for each temperature and wavenumber as
/// the integral across a wavenumber interval, as a spectral
/// irradiance in W m-2
void
planck_function(const adept::Vector& temperature,       ///< Temperature in K
		const adept::Vector& wavenumber_cm_1,   ///< Wavenumber in cm-1
		const adept::Vector& d_wavenumber_cm_1, ///< Wavenumber interval in cm-1
		adept::Matrix ans                       ///< Planck function in W m-2
		);

/// Compute the Planck function for a single temperature
inline
void
planck_function(adept::Real temperature,                ///< Temperature in K
		const adept::Vector& wavenumber_cm_1,   ///< Wavenumber in cm-1
		const adept::Vector& d_wavenumber_cm_1, ///< Wavenumber interval in cm-1
		adept::Vector ans                       ///< Planck function in W m-2
		) {
  using namespace adept;
  adept::Matrix tmp(1,wavenumber_cm_1.size());
  adept::Vector temp(1);
  temp(0) = temperature;
  planck_function(temp, wavenumber_cm_1, d_wavenumber_cm_1, tmp);
  ans = tmp(0,__);  
}
