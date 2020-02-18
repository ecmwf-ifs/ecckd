#ifndef HEATING_RATE_H
#define HEATING_RATE_H 1

#include <adept_arrays.h>

#include "constants.h"

/// Compute spectral heating rate profile from profiles of upwelling
/// and downwelling flux profiles, where fluxes and heating rates are
/// dimensioned (pressure, wavenumber)
template <bool IsActive>
void
heating_rate(const adept::Vector& pressure,                      ///< Half-level pressure (Pa)
	     const adept::Array<2,adept::Real,IsActive> flux_dn, ///< Spectral flux down (W m-2)
	     const adept::Array<2,adept::Real,IsActive> flux_up, ///< Spectral flux up (W m-2)
	     adept::Array<2,adept::Real,IsActive> hr             ///< Spectral heating rate (K s-1)
	     ) {

  using namespace adept;

  int nwav = flux_dn.size(1);

  // Factor to convert from difference in net flux across a layer to heating rate
  Vector conversion = -(ACCEL_GRAVITY/SPECIFIC_HEAT_AIR) / (pressure(range(1,end))-pressure(range(0,end-1)));
  hr = spread<1>(conversion,nwav) * (  flux_dn(range(1,end),__) - flux_dn(range(0,end-1),__)
				       - flux_up(range(1,end),__) + flux_up(range(0,end-1),__)  );
}

/// As above but for a single profile
template <bool IsActive>
void
heating_rate_single(const adept::Vector& pressure,               ///< Half-level pressure (Pa)
	     const adept::Array<1,adept::Real,IsActive> flux_dn, ///< Flux down (W m-2)
	     const adept::Array<1,adept::Real,IsActive> flux_up, ///< Flux up (W m-2)
	     adept::Array<1,adept::Real,IsActive> hr             ///< Heating rate (K s-1)
	     ) {

  using namespace adept;

  // Factor to convert from difference in net flux across a layer to heating rate
  hr = -((ACCEL_GRAVITY/SPECIFIC_HEAT_AIR) / (pressure(range(1,end))-pressure(range(0,end-1))))
    * (  flux_dn(range(1,end)) - flux_dn(range(0,end-1))
	 - flux_up(range(1,end)) + flux_up(range(0,end-1))  );
}

#endif
