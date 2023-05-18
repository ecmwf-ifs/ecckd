// average_optical_depth.cpp - Average optical depth in a spectral interval
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
				 Matrix molar_abs,            ///< Molar absorption coefficient, m2 mol-1 (pressure,g-point)
				 Matrix min_molar_abs,        ///< Minimum molar absorption coefficient
				 Matrix max_molar_abs)        ///< Maximum molar absorption coefficient

{
  const Real OD_SCALING = 1.0;

#pragma omp parallel for
  for (int ig = 0; ig < ng; ++ig) {
    // Fails if index is empty, due to earlier problems...
    intVector index = find(g_point == ig);
    Vector optical_depth_fit, min_optical_depth, max_optical_depth;

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
    else if (averaging_method == "transmission-3") {
      optical_depth_fit = min(0.9999999999999999,sum((1.0-exp(-optical_depth.soft_link()(__,index)*(LW_DIFFUSIVITY*OD_SCALING*3.0)))
					    * planck_fl.soft_link()(__,index), 1)
			    / sum(planck_fl.soft_link()(__,index),1));
      optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*OD_SCALING*3.0));
    }
    else if (averaging_method == "transmission-10") {
      optical_depth_fit = min(0.9999999999999999,sum((1.0-exp(-optical_depth.soft_link()(__,index)*(LW_DIFFUSIVITY*OD_SCALING*10.0)))
					    * planck_fl.soft_link()(__,index), 1)
			    / sum(planck_fl.soft_link()(__,index),1));
      optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*OD_SCALING*10.0));
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
    else if (averaging_method == "hybrid-logarithmic-transmission-3") {
      optical_depth_fit.resize(optical_depth.size(0));
      for (int iz = 0; iz < optical_depth.size(0); ++iz) {
	if (pressure_fl(iz) > 100.0e2) {
	  // Logarithmic averaging for pressure larger than 100 hPa
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
	else {
	  // Transmission-3 averaging for smaller pressures
	  optical_depth_fit(iz) = std::min(0.9999999999999999,sum((1.0-exp(-optical_depth.soft_link()(iz,index)*(LW_DIFFUSIVITY*OD_SCALING*3.0)))
							     * planck_fl.soft_link()(iz,index))
				      / sum(planck_fl.soft_link()(iz,index)));
	  optical_depth_fit(iz) = abs(-log(1.0-optical_depth_fit(iz))/(LW_DIFFUSIVITY*OD_SCALING*3.0));
	}
      }
    }
    else {
      ERROR << "averaging_method \"" << averaging_method << "\" not understood";
      THROW(PARAMETER_ERROR);
    }

    if (index.empty()) {
      WARNING << "No wavenumbers with g_point == " << ig << ": skipping\n";
      ENDWARNING;
      optical_depth_fit = 0.0;
      min_optical_depth = 0.0;
      max_optical_depth = 0.0;
    }
    else {
      // The calculations above can lead to the fitted optical depth
      // not lying between min and max, which needs to be fixed or the
      // bounded optimization will fail
      min_optical_depth = minval(optical_depth.soft_link()(__,index), 1);
      max_optical_depth = maxval(optical_depth.soft_link()(__,index), 1);

      // Sometimes transmission averaging goes inaccurate for low
      // optical depths: ensure it is in the range max to min
      optical_depth_fit = max(min_optical_depth, min(optical_depth_fit, max_optical_depth));
      
      if (any(min_optical_depth > optical_depth_fit)) {
	WARNING << "G-point " << ig << ": min optical depth > average optical depth: correcting; length(index) = " << index.size() << "\n";
	ENDWARNING;
	min_optical_depth.where(min_optical_depth > optical_depth_fit) = optical_depth_fit;
      }
      if (any(min_optical_depth > 0.0 && min_optical_depth >= max_optical_depth)) {
	WARNING << "G-point " << ig << ": min optical depth >= max_optical_depth: length(index) = " << index.size() << "\n";
	ENDWARNING;
	// Bounded minimization can't cope with equal bounds
	intVector index2 = find(min_optical_depth > 0.0 && min_optical_depth >= max_optical_depth);
	min_optical_depth(index2) *= 0.99;
	max_optical_depth(index2) *= 1.01;
      }
    }

    // Abs needed because -log(1) is -0 and we want to remove the sign
    // from a negative zero.
    if (reference_surface_vmr > 0.0) {
      // Molar absorption coefficient, where reference_surface_vmr
      // actually refers to the volume mixing ratio at all heights in
      // this Idealized dataset
      molar_abs.soft_link()(__,ig) 
	= ((ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR) / reference_surface_vmr)
	* optical_depth_fit / (pressure_hl.soft_link()(range(1,end))-pressure_hl.soft_link()(range(0,end-1)));
      if (!min_molar_abs.empty()) {
	min_molar_abs.soft_link()(__,ig)
	  = ((ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR) / reference_surface_vmr)
	  * min_optical_depth / (pressure_hl.soft_link()(range(1,end))-pressure_hl.soft_link()(range(0,end-1)));
	max_molar_abs.soft_link()(__,ig)
	  = ((ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR) / reference_surface_vmr)
	  * max_optical_depth / (pressure_hl.soft_link()(range(1,end))-pressure_hl.soft_link()(range(0,end-1)));
      }
    }
    else {
      // Calculate mean optical depth, rather than molar absorption
      // coefficient
      molar_abs.soft_link()(__,ig)     = optical_depth_fit;
      if (!min_molar_abs.empty()) {
	min_molar_abs.soft_link()(__,ig) = min_optical_depth;
	max_molar_abs.soft_link()(__,ig) = max_optical_depth;
      }
    }
  }

}
