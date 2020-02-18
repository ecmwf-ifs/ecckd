#include "constants.h"

using namespace adept;

//#define LINEAR 1

// Average optical depths for each pressure level to each g point
void
average_optical_depth_to_g_point(int ng,                      ///< Number of g points
				 Real reference_surface_vmr,  ///< Volume mixing ratio (mol mol-1)
				 const Vector& pressure_fl,   ///< Full-level pressure (Pa)
				 const Vector& pressure_hl,   ///< Half-level pressure (Pa)
				 const Vector& g_point,       ///< G point for each wavenumber
				 const Matrix& optical_depth, ///< Optical depth (pressure,wavenumber)
				 const Matrix& planck_fl,     ///< Planck function, W m-2 (pressure,wavenumber)
				 Matrix molar_abs)            ///< Molar absorption coefficient, m2 mol-1 (pressure,g-point)
{
  const Real OD_SCALING = 1.0;

#pragma omp parallel for
  for (int ig = 0; ig < ng; ++ig) {
    // Fails if index is empty, due to earlier problems...
    intVector index = find(g_point == ig);
    Vector optical_depth_fit;

#ifdef LINEAR
    optical_depth_fit = sum(optical_depth.soft_link()(__,index)*planck_fl.soft_link()(__,index), 1)
      / sum(planck_fl.soft_link()(__,index),1);
#else
    optical_depth_fit = min(0.9999999999999999,sum((1.0-exp(-optical_depth.soft_link()(__,index)*(LW_DIFFUSIVITY*OD_SCALING)))
					    * planck_fl.soft_link()(__,index), 1)
			    / sum(planck_fl.soft_link()(__,index),1));
    optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*OD_SCALING));
#endif
    // Abs needed because -log(1) is -0 and we want to remove the sign
    // from a negative zero.
    if (reference_surface_vmr > 0.0) {
      molar_abs.soft_link()(__,ig) 
	= ((ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR) / reference_surface_vmr)
	* optical_depth_fit / (pressure_hl.soft_link()(range(1,end))-pressure_hl.soft_link()(range(0,end-1)));
    }
    else {
      // Calculate mean optical depth, rather than molar absorption
      // coefficient
      molar_abs.soft_link()(__,ig) = optical_depth_fit;
    }
  }

}
