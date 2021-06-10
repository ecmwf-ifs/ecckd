// radiative_transfer_sw.cpp - Perform shortwave radiative transfer
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

#include "radiative_transfer_sw.h"

/// Perform direct-beam shortwave radiative transfer using solar
/// irradiance integrated across each spectral interval and optical
/// depth (dimensioned layer,wavenumber) to obtain upwelling and
/// downwelling spectral fluxes (dimensioned
/// half-level,wavenumber). The function is templated so that
/// optionally automatic differentiation can be used to compute the
/// derivative of the fluxes with respect to the optical depth.
template <bool IsActive>
void
radiative_transfer_direct_sw(adept::Real cos_sza,             ///< Cosine of the solar zenith angle
			     const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
			     const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
			     adept::Array<2,adept::Real,IsActive> flux_dn ///< Spectral flux down in W m-2
			     ) {
  using namespace adept;

  int nlay = optical_depth.size(0);
  //int nwav = optical_depth.size(1);

  Real minus_sec_sza = -1.0 / cos_sza;

  // Work down from top of atmosphere
  flux_dn(0,__) = cos_sza*ssi;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux_dn(ilay+1,__) = flux_dn(ilay,__) * exp(minus_sec_sza*optical_depth(ilay,__));
  }
}

/// As radiative_transfer_direct_sw but the upwelling flux is also
/// computed assuming no Rayleigh scattering
template <bool IsActive>
void
radiative_transfer_norayleigh_sw(adept::Real cos_sza,      ///< Cosine of the solar zenith angle
				 const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
				 const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
				 const adept::Vector& albedo, ///< Surface albedo
				 adept::Array<2,adept::Real,IsActive> flux_dn, ///< Spectral flux down in W m-2
				 adept::Array<2,adept::Real,IsActive> flux_up ///< Spectral flux up in W m-2
				 ) {
  using namespace adept;

  int nlay = optical_depth.size(0);
  //int nwav = optical_depth.size(1);

  Real minus_sec_sza = -1.0 / cos_sza;

  // Negative of the secant of the two-stream zenith angle, which in
  // the shortwave we take as 60 degrees consistent with Zdunkowski
  // (1980)
  static const Real minus_sec_tsza = -2.0;

  // Work down from top of atmosphere
  flux_dn(0,__) = cos_sza*ssi;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux_dn(ilay+1,__) = flux_dn(ilay,__) * exp(minus_sec_sza*optical_depth(ilay,__));
  }
  flux_up(nlay,__) = flux_dn(nlay,__) * albedo;
  for (int ilay = nlay-1; ilay >= 0; --ilay) {
    flux_up(ilay,__) = flux_up(ilay+1,__) * exp(minus_sec_tsza*optical_depth(ilay,__));
  }
}

// Explicit instantiations
template
void
radiative_transfer_direct_sw<true>(adept::Real cos_sza,             ///< Cosine of the solar zenith angle
			    const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
			    const adept::Array<2,adept::Real,true>& optical_depth, ///< Layer optical depth
			    adept::Array<2,adept::Real,true> flux_dn ///< Spectral flux down in W m-2
			    );

template
void
radiative_transfer_direct_sw<false>(adept::Real cos_sza,             ///< Cosine of the solar zenith angle
			     const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
			     const adept::Array<2,adept::Real,false>& optical_depth, ///< Layer optical depth
			     adept::Array<2,adept::Real,false> flux_dn ///< Spectral flux down in W m-2
			     );
template
void
radiative_transfer_norayleigh_sw<true>(adept::Real cos_sza,      ///< Cosine of the solar zenith angle
				 const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
				 const adept::Array<2,adept::Real,true>& optical_depth, ///< Layer optical depth
				 const adept::Vector& albedo, ///< Surface albedo
				 adept::Array<2,adept::Real,true> flux_dn, ///< Spectral flux down in W m-2
				 adept::Array<2,adept::Real,true> flux_up ///< Spectral flux up in W m-2
				       );
template
void
radiative_transfer_norayleigh_sw<false>(adept::Real cos_sza,      ///< Cosine of the solar zenith angle
				 const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
				 const adept::Array<2,adept::Real,false>& optical_depth, ///< Layer optical depth
				 const adept::Vector& albedo, ///< Surface albedo
				 adept::Array<2,adept::Real,false> flux_dn, ///< Spectral flux down in W m-2
				 adept::Array<2,adept::Real,false> flux_up ///< Spectral flux up in W m-2
				       );

/// As radiative_transfer_direct_sw but with broadband flux outputs, to
/// reduce memory consumption if not needed
template <bool IsActive>
void
radiative_transfer_direct_sw_bb(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
			 const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
			 const adept::Array<2,adept::Real,IsActive>& spectral_od, ///< Spectral layer optical depth
			 const adept::Array<1,adept::Real,IsActive>& grey_od,     ///< Grey layer optical depth
			 adept::Array<1,adept::Real,IsActive> flux_dn ///< Broadband flux down in W m-2
			 ) {
  using namespace adept;

  int nlay = spectral_od.size(0);
  int nwav = spectral_od.size(1);

  Real minus_sec_sza = -1.0 / cos_sza;

  // Spectral downward flux
  adept::Array<1,Real,IsActive> flux(nwav);

  // Work down from top of atmosphere: here "flux" is the spectral downward flux
  flux_dn(0) = cos_sza*sum(ssi);
  flux = cos_sza*ssi;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux = flux * exp(minus_sec_sza*(spectral_od(ilay,__) + grey_od(ilay)));
    flux_dn(ilay+1) = sum(flux);
  }
}

/// As radiative_transfer_norayleigh_sw but with broadband flux
/// outputs, to reduce memory consumption if not needed
template <bool IsActive>
void
radiative_transfer_norayleigh_sw_bb(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
				    const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
				    const adept::Array<2,adept::Real,IsActive>& spectral_od, ///< Spectral layer optical depth
				    const adept::Array<1,adept::Real,IsActive>& grey_od,     ///< Grey layer optical depth
				    adept::Real albedo, ///< Surface albedo
				    adept::Array<1,adept::Real,IsActive> flux_dn, ///< Broadband flux down in W m-2
				    adept::Array<1,adept::Real,IsActive> flux_up ///< Broadband flux up in W m-2
			 ) {
  using namespace adept;

  int nlay = spectral_od.size(0);
  int nwav = spectral_od.size(1);

  Real minus_sec_sza = -1.0 / cos_sza;

  // Negative of the secant of the two-stream zenith angle, which in
  // the shortwave we take as 60 degrees consistent with Zdunkowski
  // (1980)
  static const Real minus_sec_tsza = -2.0;

  // Spectral downward/upward flux
  adept::Array<1,Real,IsActive> flux(nwav);

  // Work down from top of atmosphere: here "flux" is the spectral downward flux
  flux_dn(0) = cos_sza*sum(ssi);
  flux = cos_sza*ssi;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux = flux * exp(minus_sec_sza*(spectral_od(ilay,__) + grey_od(ilay)));
    flux_dn(ilay+1) = sum(flux);
  }
  // Now flux is upward
  flux *= albedo;
  flux_up(nlay) = sum(flux);
  for (int ilay = nlay-1; ilay >= 0; --ilay) {
    flux = flux * exp(minus_sec_tsza*(spectral_od(ilay,__) + grey_od(ilay)));
    flux_up(ilay) = sum(flux);
  }
}

// Explicit instantiations
template
void
radiative_transfer_direct_sw_bb<true>(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
			 const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
			 const adept::Array<2,adept::Real,true>& spectral_od, ///< Spectral layer optical depth
			 const adept::Array<1,adept::Real,true>& grey_od,     ///< Grey layer optical depth
			 adept::Array<1,adept::Real,true> flux_dn ///< Broadband flux down in W m-2
			 );
template
void
radiative_transfer_direct_sw_bb<false>(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
			 const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
			 const adept::Array<2,adept::Real,false>& spectral_od, ///< Spectral layer optical depth
			 const adept::Array<1,adept::Real,false>& grey_od,     ///< Grey layer optical depth
			 adept::Array<1,adept::Real,false> flux_dn ///< Broadband flux down in W m-2
			 );

template
void
radiative_transfer_norayleigh_sw_bb<true>(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
			 const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
			 const adept::Array<2,adept::Real,true>& spectral_od, ///< Spectral layer optical depth
			 const adept::Array<1,adept::Real,true>& grey_od,     ///< Grey layer optical depth
			 adept::Real albedo, ///< Surface albedo
					  adept::Array<1,adept::Real,true> flux_dn, ///< Broadband flux down in W m-2
			 adept::Array<1,adept::Real,true> flux_up ///< Broadband flux up in W m-2
			 );

template
void
radiative_transfer_norayleigh_sw_bb<false>(adept::Real cos_sza,         ///< Cosine of the solar zenith angle
			 const adept::Vector& ssi,    ///< Spectral solar irradiance in W m-2
			 const adept::Array<2,adept::Real,false>& spectral_od, ///< Spectral layer optical depth
			 const adept::Array<1,adept::Real,false>& grey_od,     ///< Grey layer optical depth
			 adept::Real albedo, ///< Surface albedo
			 adept::Array<1,adept::Real,false> flux_dn, ///< Broadband flux down in W m-2
			 adept::Array<1,adept::Real,false> flux_up ///< Broadband flux up in W m-2
			 );
