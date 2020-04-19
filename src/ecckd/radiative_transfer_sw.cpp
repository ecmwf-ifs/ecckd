#include "radiative_transfer_sw.h"

/// Perform direct-beam shortwave radiative transfer using planck function
/// integrated across each spectral interval (dimensioned
/// half-level,wavenumber) and optical depth (dimensioned
/// layer,wavenumber) to obtain upwelling and downwelling spectral
/// fluxes (dimensioned half-level,wavenumber). The function is
/// templated so that optionally automatic differentiation can be used
/// to compute the derivative of the fluxes with respect to the
/// optical depth.
template <bool IsActive>
void
radiative_transfer_direct_sw(adept::Real cos_sza,             ///< Cosine of the solar zenith angle
			     const adept::Vector& ssi, ///< Spectral solar irradiance in W m-2
			     const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
			     adept::Array<2,adept::Real,IsActive> flux_dn ///< Spectral flux down in W m-2
			     ) {
  using namespace adept;

  int nlay = optical_depth.size(0);
  int nwav = optical_depth.size(1);

  Real minus_sec_sza = -1.0 / cos_sza;

  // Work down from top of atmosphere
  flux_dn(0,__) = ssi;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux_dn(ilay+1,__) = flux_dn(ilay,__) * exp(minus_sec_sza*optical_depth(ilay,__));
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

/// As radiative_transfer_sw but with broadband flux outputs, to
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
  flux_dn(0) = sum(ssi);
  flux = ssi;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux = flux * exp(minus_sec_sza*(spectral_od(ilay,__) + grey_od(ilay)));
    flux_dn(ilay+1) = sum(flux);
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
