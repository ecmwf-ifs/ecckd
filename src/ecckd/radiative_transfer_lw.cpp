#include "radiative_transfer_lw.h"

/// Perform longwave radiative transfer using planck function
/// integrated across each spectral interval (dimensioned
/// half-level,wavenumber) and optical depth (dimensioned
/// layer,wavenumber) to obtain upwelling and downwelling spectral
/// fluxes (dimensioned half-level,wavenumber). The function is
/// templated so that optionally automatic differentiation can be used
/// to compute the derivative of the fluxes with respect to the
/// optical depth.
template <bool IsActive>
void
radiative_transfer_lw(const adept::Matrix& planck,                  ///< Planck function in W m-2
		      const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
		      const adept::Vector& surf_emissivity,         ///< Surface emissivity
		      const adept::Vector& surf_planck,             ///< Surface Planck function in W m-2
		      adept::Array<2,adept::Real,IsActive> flux_dn, ///< Spectral flux down in W m-2
		      adept::Array<2,adept::Real,IsActive> flux_up  ///< Spectral flux up in W m-2
		      ) {
  using namespace adept;

  int nlay = optical_depth.size(0);
  int nwav = optical_depth.size(1);

  adept::Array<2,Real,IsActive> emissivity(nlay,nwav), factor(nlay,nwav);

  emissivity = 1.0 - exp(-LW_DIFFUSIVITY*optical_depth);
  factor.where(emissivity > 1.0e-5) = either_or(1.0 - emissivity*(1.0/LW_DIFFUSIVITY)/optical_depth,
						0.5 * emissivity);
  // Work down from top of atmosphere
  flux_dn(0,__) = 0.0;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    flux_dn(ilay+1,__) = flux_dn(ilay,__) * (1.0 - emissivity(ilay,__))
      + planck(ilay,__)   * (emissivity(ilay,__)-factor(ilay,__))
      + planck(ilay+1,__) * factor(ilay,__);
  }

  flux_up(nlay,__) = surf_planck * surf_emissivity
    + (1.0 - surf_emissivity) * flux_dn(nlay,__);
  // Work up from surface
  for (int ilay = nlay-1; ilay >= 0; --ilay) {
    flux_up(ilay,__) = flux_up(ilay+1,__) * (1.0 - emissivity(ilay,__))
      + planck(ilay+1,__) * (emissivity(ilay,__)-factor(ilay,__))
      + planck(ilay,__)   * factor(ilay,__);
  }
}

// Explicit instantiations
template
void
radiative_transfer_lw<true>(const adept::Matrix& planck,              ///< Planck function in W m-2
			    const adept::Array<2,adept::Real,true>& optical_depth, ///< Layer optical depth
			    const adept::Vector& surf_emissivity,     ///< Surface emissivity
			    const adept::Vector& surf_planck,         ///< Surface Planck function in W m-2
			    adept::Array<2,adept::Real,true> flux_dn, ///< Spectral flux down in W m-2
			    adept::Array<2,adept::Real,true> flux_up  ///< Spectral flux up in W m-2
			    );
template
void
radiative_transfer_lw<false>(const adept::Matrix& planck,               ///< Planck function in W m-2
			     const adept::Array<2,adept::Real,false>& optical_depth, ///< Layer optical depth
			     const adept::Vector& surf_emissivity,      ///< Surface emissivity
			     const adept::Vector& surf_planck,          ///< Surface Planck function in W m-2
			     adept::Array<2,adept::Real,false> flux_dn, ///< Spectral flux down in W m-2
			     adept::Array<2,adept::Real,false> flux_up  ///< Spectral flux up in W m-2
			     );


/// As radiative_transfer_lw but with broadband flux outputs, to
/// reduce memory consumption if not needed
template <bool IsActive>
void
radiative_transfer_lw_bb(const adept::Matrix& planck,                  ///< Planck function in W m-2
			 const adept::Array<2,adept::Real,IsActive>& optical_depth, ///< Layer optical depth
			 const adept::Vector& surf_emissivity,         ///< Surface emissivity
			 const adept::Vector& surf_planck,             ///< Surface Planck function in W m-2
			 adept::Array<1,adept::Real,IsActive> flux_dn, ///< Broadband flux down in W m-2
			 adept::Array<1,adept::Real,IsActive> flux_up  ///< Broadband flux up in W m-2
			 ) {
  using namespace adept;

  int nlay = optical_depth.size(0);
  int nwav = optical_depth.size(1);

  // 1D only to conserve memory
  adept::Array<1,Real,IsActive> emissivity(nwav), factor(nwav);

  // Spectral downward or upward flux
  adept::Array<1,Real,IsActive> flux(nwav);

  // Work down from top of atmosphere: here "flux" is the spectral downward flux
  flux = 0.0;
  flux_dn(0) = 0.0;
  for (int ilay = 0; ilay < nlay; ++ilay) {
    emissivity = 1.0 - exp(-LW_DIFFUSIVITY*optical_depth(ilay,__));
    factor.where(emissivity > 1.0e-5) = either_or(1.0 - emissivity*(1.0/LW_DIFFUSIVITY)/optical_depth(ilay,__),
						  0.5 * emissivity);
    flux = flux * (1.0 - emissivity)
      + planck(ilay,__)   * (emissivity-factor)
      + planck(ilay+1,__) * factor;
    flux_dn(ilay+1) = sum(flux);
  }

  // Flux is now the spectral upward flux
  flux = surf_planck * surf_emissivity + (1.0 - surf_emissivity) * flux;
  flux_up(nlay) = sum(flux);
  // Work up from surface
  for (int ilay = nlay-1; ilay >= 0; --ilay) {
    emissivity = 1.0 - exp(-LW_DIFFUSIVITY*optical_depth(ilay,__));
    factor.where(emissivity > 1.0e-5) = either_or(1.0 - emissivity*(1.0/LW_DIFFUSIVITY)/optical_depth(ilay,__),
						  0.5 * emissivity);
    flux = flux * (1.0 - emissivity)
      + planck(ilay+1,__) * (emissivity-factor)
      + planck(ilay,__)   * factor;
    flux_up(ilay) = sum(flux);
  }
}

// Explicit instantiations
template
void
radiative_transfer_lw_bb<true>(const adept::Matrix& planck,              ///< Planck function in W m-2
			    const adept::Array<2,adept::Real,true>& optical_depth, ///< Layer optical depth
			    const adept::Vector& surf_emissivity,     ///< Surface emissivity
			    const adept::Vector& surf_planck,         ///< Surface Planck function in W m-2
			    adept::Array<1,adept::Real,true> flux_dn, ///< Broadband flux down in W m-2
			    adept::Array<1,adept::Real,true> flux_up  ///< Broadband flux up in W m-2
			    );
template
void
radiative_transfer_lw_bb<false>(const adept::Matrix& planck,               ///< Planck function in W m-2
			     const adept::Array<2,adept::Real,false>& optical_depth, ///< Layer optical depth
			     const adept::Vector& surf_emissivity,      ///< Surface emissivity
			     const adept::Vector& surf_planck,          ///< Surface Planck function in W m-2
			     adept::Array<1,adept::Real,false> flux_dn, ///< Broadband flux down in W m-2
			     adept::Array<1,adept::Real,false> flux_up  ///< Broadband flux up in W m-2
			     );
