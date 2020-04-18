#ifndef AVERAGE_OPTICAL_DEPTH_H
#define AVERAGE_OPTICAL_DEPTH_H

// Average optical depths for each pressure level to each g point
void
average_optical_depth_to_g_point(int ng,                             ///< Number of g points
				 adept::Real reference_surface_vmr,  ///< Volume mixing ratio (mol mol-1)
				 const adept::Vector& pressure_fl,   ///< Full-level pressure (Pa)
				 const adept::Vector& pressure_hl,   ///< Half-level pressure (Pa)
				 const adept::Vector& g_point,       ///< G point for each wavenumber
				 const adept::Matrix& optical_depth, ///< Optical depth (pressure,wavenumber)
				 const adept::Matrix& planck_fl,     ///< Planck function, W m-2 (pressure,wavenumber)
				 const std::string& averaging_method,
				 adept::Matrix molar_abs);           ///< Molar absorption coefficient, m2 mol-1 (pressure,g-point)
#endif
