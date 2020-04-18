#include "constants.h"
#include "Logging.h"

using namespace adept;

// Average optical depths for each pressure level to each g point
void
average_optical_depth_to_g_point(int ng,                      ///< Number of g points
				 Real reference_surface_vmr,  ///< Volume mixing ratio (mol mol-1)
				 const Vector& pressure_fl,   ///< Full-level pressure (Pa)
				 const Vector& pressure_hl,   ///< Half-level pressure (Pa)
				 const Vector& g_point,       ///< G point for each wavenumber
				 const Matrix& optical_depth, ///< Optical depth (pressure,wavenumber)
				 const Matrix& planck_fl,     ///< Planck function, W m-2 (pressure,wavenumber)
				 const std::string& averaging_method,
				 Matrix molar_abs)            ///< Molar absorption coefficient, m2 mol-1 (pressure,g-point)

{
  const Real OD_SCALING = 1.0;

#pragma omp parallel for
  for (int ig = 0; ig < ng; ++ig) {
    // Fails if index is empty, due to earlier problems...
    intVector index = find(g_point == ig);
    Vector optical_depth_fit;

    if (averaging_method == "linear") {
      optical_depth_fit = sum(optical_depth.soft_link()(__,index)*planck_fl.soft_link()(__,index), 1)
	/ sum(planck_fl.soft_link()(__,index),1);
    }
    else if (averaging_method == "transmission") {
      optical_depth_fit = min(0.9999999999999999,sum((1.0-exp(-optical_depth.soft_link()(__,index)*(LW_DIFFUSIVITY*OD_SCALING)))
					    * planck_fl.soft_link()(__,index), 1)
			    / sum(planck_fl.soft_link()(__,index),1));
      optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*OD_SCALING));
    }
    else if (averaging_method == "transmission-2") {
      optical_depth_fit = min(0.9999999999999999,sum((1.0-exp(-optical_depth.soft_link()(__,index)*(LW_DIFFUSIVITY*OD_SCALING*2.0)))
					    * planck_fl.soft_link()(__,index), 1)
			    / sum(planck_fl.soft_link()(__,index),1));
      optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*OD_SCALING*2.0));
    }
    else if (averaging_method == "square-root") {
      optical_depth_fit = sum(sqrt(optical_depth.soft_link()(__,index))*planck_fl.soft_link()(__,index), 1)
	/ sum(planck_fl.soft_link()(__,index),1);
      optical_depth_fit *= optical_depth_fit;
    }
    else if (averaging_method == "logarithmic") {
      optical_depth_fit.resize(optical_depth.size(0));
      for (int iz = 0; iz < optical_depth.size(0); ++iz) {
	intVector iindex = find(g_point == ig && optical_depth.soft_link()(iz,__) > 0.0);
        Vector od_nonzero = optical_depth.soft_link()(iz,iindex);
        if (od_nonzero.size() == index.size()) {
          // Pure logarithmic average
          optical_depth_fit(iz) = exp(sum(log(od_nonzero)*planck_fl.soft_link()(iz,index))
				      / sum(planck_fl.soft_link()(iz,index)));
        }
	else if (od_nonzero.empty()) {
	  // No non-zero data
	  optical_depth_fit(iz) = 0.0;
	}
	else {
	  // Some zeros: logarithmic average of non-zeros then linear average with zeros
          optical_depth_fit(iz) = exp(sum(log(od_nonzero)*planck_fl.soft_link()(iz,iindex))
				      / sum(planck_fl.soft_link()(iz,iindex)))
	    * (static_cast<Real>(iindex.size())/static_cast<Real>(index.size()));
	}
      }
    }
    else {
      ERROR << "averaging_method \"" << averaging_method << "\" not understood";
      THROW(PARAMETER_ERROR);
    }
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
