// ckd_model.h - Class for storing a CKD model
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

#ifndef CKD_MODEL_H
#define CKD_MODEL_H 1

#include <vector>
#include <string>

#include <adept_arrays.h>

#include "Error.h"
#include "rayleigh_scattering.h"
#include "calc_cost_function_sw.h"

using namespace adept;

// 0=none, 1=linear, 2=look-up table
typedef enum {
  NONE = 0,
  LINEAR,
  LUT,
  RELATIVE_LINEAR,
} ConcDependence;

template <bool IsActive>
struct SingleGasData {
  SingleGasData(const std::string& molecule_)
    : molecule(molecule_) {
    Molecule = molecule;
    std::transform(Molecule.begin(), Molecule.end(), Molecule.begin(), ::toupper);
  }

  std::string molecule, Molecule;

  ConcDependence conc_dependence;

  // Molar absorption coefficient in m2 mol-1. If conc_dependence=0
  // then it is the absorption cross section per mole of dry air.  If
  // conc_dependence=1, it is the absorption cross section per mole of
  // the gas in question. It is dimensioned
  // (temperature,pressure,g_point).
  Array<3,adept::Real,IsActive> molar_abs;

  // If conc_dependence=2, then we have an additional dimension for
  // concentration.  It is dimensioned
  // (conc,temperature,pressure,g_point).
  Array<4,adept::Real,IsActive> molar_abs_conc;

  // If conc_dependence=3 then the following reference concentration
  // is subtracted from the actual concentration before the result is
  // multiplied by the mass aborption coefficient
  Real reference_vmr;

  // Volume mixing ratio coordinate variable
  Vector vmr;

  // Inverse of the background error covariance matrix for one g point
  SparseMatrix inv_background;
  //SymmMatrix inv_background;

  // Index of the first element of the full state variable x in the
  // present gas
  int ix;

  // Is this gas to be optimized?
  bool is_active = false;

  // Information about single gases contributing to a composite gas
  Matrix composite_vmr;
  std::string composite_molecules;

};

/// Data describing a correlated k-distribution model
template <bool IsActive>
class CkdModel {

public:

  /// Construct a CKD model from a file
  CkdModel(const std::string& file_name,
	   const std::vector<std::string>& gas_list = std::vector<std::string>()) {
    read(file_name, gas_list);
  }

  /// Longwave constructor
  CkdModel(const std::vector<SingleGasData<false> >& single_gas_data,
	   const Vector& temperature_planck,
	   const Matrix& planck_function,
	   const Vector& pressure, 
	   const Matrix& temperature,
	   const Vector& wavenumber1,
	   const Vector& wavenumber2,
	   const Matrix& gpoint_fraction,
	   const Vector& wavenumber1_band,
	   const Vector& wavenumber2_band,
	   const intVector& band_number,
	   const std::string& history = std::string(),
	   const std::string& config = std::string()) 
    : single_gas_data_(single_gas_data),
      temperature_planck_(eval(temperature_planck)),
      planck_function_(eval(planck_function)),
      log_pressure_(log(pressure)),
      temperature_(eval(temperature)),
      wavenumber1_(eval(wavenumber1)),
      wavenumber2_(eval(wavenumber2)),
      gpoint_fraction_(eval(gpoint_fraction)),
      wavenumber1_band_(eval(wavenumber1_band)),
      wavenumber2_band_(eval(wavenumber2_band)),
      band_number_(eval(band_number)),
      history_(history),
      config_(config) {
    ng_ = planck_function.size(1);
    nt_ = temperature.size(0);
    np_ = pressure.size();
    nwav_ = wavenumber1.size();
  }

  /// Shortwave constructor
  CkdModel(const std::vector<SingleGasData<false> >& single_gas_data,
	   const Vector& solar_irradiance, // per g point
	   const Vector& pressure, 
	   const Matrix& temperature,
	   const Vector& wavenumber1,
	   const Vector& wavenumber2,
	   const Matrix& gpoint_fraction,
	   const Vector& ssi_intervals, // In the intervals bounded by wavenumber1,2
	   const Vector& wavenumber1_band,
	   const Vector& wavenumber2_band,
	   const intVector& band_number,
	   const std::string& history = std::string(),
	   const std::string& config = std::string()) 
    : single_gas_data_(single_gas_data),
      solar_irradiance_(eval(solar_irradiance)),
      log_pressure_(log(pressure)),
      temperature_(eval(temperature)),
      wavenumber1_(eval(wavenumber1)),
      wavenumber2_(eval(wavenumber2)),
      gpoint_fraction_(eval(gpoint_fraction)),
      wavenumber1_band_(eval(wavenumber1_band)),
      wavenumber2_band_(eval(wavenumber2_band)),
      band_number_(eval(band_number)),
      history_(history),
      config_(config) {
    ng_ = solar_irradiance.size();
    nt_ = temperature.size(0);
    np_ = pressure.size();
    nwav_ = wavenumber1.size();
    calc_rayleigh_molar_scat(ssi_intervals);
  }

  /// Read a CKD model from a file
  void read(const std::string& file_name,
	    const std::vector<std::string>& gas_list = std::vector<std::string>());

  // Write a CKD model to file
  void write(const std::string& file_name, int argc, const char* argv[], const std::string& config_str);

  /// Calculate the optical depth at each g-point for multiple
  /// atmospheric columns, where the arguments are dimensioned
  /// (column,level) and the output dimensioned (column,level,g-point)
  Array<3,Real,IsActive> calc_optical_depth(int igas,         ///< Gas number 
			    const Matrix& pressure_hl,        ///< Pressure at half levels (Pa)
			    const Matrix& temperature_fl,     ///< Temperature at full levels (K)
			    const Matrix& vmr_fl = Matrix()); ///< Volume mixing ratio at full levels

  /// As above but provide the name of the gas
  Array<3,Real,IsActive> calc_optical_depth(const std::string& gas, ///< Gas name, lower case
			    const Matrix& pressure_hl,        ///< Pressure at half levels (Pa)
			    const Matrix& temperature_fl,     ///< Temperature at full levels (K)
			    const Matrix& vmr_fl = Matrix()) ///< Volume mixing ratio at full levels
    {
      int igas = get_gas_index(gas);
      if (igas == -1) {
	ERROR << "CKD model does not contain " << gas;
	THROW(UNEXPECTED_EXCEPTION);
      }
      else {
	return calc_optical_depth(igas, pressure_hl, temperature_fl, vmr_fl);
      }
    }

  // Return the index to gas called "gas", where an empty string or
  // "composite" will match the first gas with no concentration
  // dependence
  int get_gas_index(std::string gas) {
    if (gas.empty()) {
      gas = "composite";
    }
    auto itgas = std::find(molecules.begin(), molecules.end(), gas);
    if (itgas == molecules.end()) {
      if (gas == "composite") {
	// Search for first gas 
	for (int igas = 0; igas < ngas(); ++igas) {
	  const SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
	  if (this_gas.conc_dependence == NONE) {
	    return igas;
	  }
	}
      }
      return -1; // Gas not found
    }
    else {
      return std::distance(molecules.begin(), itgas);
    }
  }

  /// Calculate the Rayleigh optical depth
  Array<3,Real,IsActive> calc_rayleigh_optical_depth(const Matrix& pressure_hl) {
    int nprof = pressure_hl.dimension(0);
    int nlev  = pressure_hl.dimension(1)-1;
    Array<3,Real,IsActive> rayleigh_od(nprof, nlev, ng_);
    Matrix moles_per_layer = (pressure_hl(__,range(1,end)) - pressure_hl(__,range(0,end-1)))
      * (1.0 / (ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR ));
    for (int ig = 0; ig < ng_; ++ig) {
      rayleigh_od(__,__,ig) = moles_per_layer * rayleigh_molar_scat_(ig);
    }
    return rayleigh_od;
  }

  /// Scale the optical depth coefficients of each gas equally, where
  /// scaling is dimensioned (nz,ng)
  void scale_optical_depth(const Vector& pressure_fl, const Matrix& scaling);

  /// Create error covariance matrices
  void create_error_covariances(Real err, Real pressure_corr, Real temperature_corr, Real conc_corr,
				Real rayleigh_prior_error = -1.0);

  /// Ensure that gases with a "relative-linear" representation cannot
  /// lead to a negative optical depth if their concentration is zero
  void cap_relative_linear_coeffts();

  /// Return the background contribution to cost function, J, and also
  /// the gradient dJ/dx, where delta_x is the difference between the
  /// current state and the prior
  Real calc_background_cost_function(const Vector& delta_x, Vector gradient);

  Array3 calc_planck_function(const Matrix& temperature_hl);
  Matrix calc_planck_function(const Vector& temperature);
  void clear() {
    single_gas_data_.clear();
    log_pressure_.clear();
    temperature_.clear();
    temperature_planck_.clear();
    planck_function_.clear();
    molecules.clear();
    x.clear();
    x_prior.clear();
  }

  /// Return list of the bands corresponding to each g-point
  intVector iband_per_g(const Vector& wavenumber1, const Vector& wavenumber2) const {
    intVector iband(ng_);
    iband = -1;
    for (int ib = 0; ib < wavenumber1.size(); ++ib) {
      Vector weight = sum(gpoint_fraction_(__,find(wavenumber1_ >= wavenumber1(ib)
						   && wavenumber2_ <= wavenumber2(ib))),1);
      if (any(weight > 0.05 && (weight < 0.95 || weight > 1.05))) {
	ERROR << "G-points do not lie entirely within requested bands: weights for band "
	      << wavenumber1(ib) << "-" << wavenumber2(ib) << " cm-1 are "
	      << weight;
	THROW(1);
      }
      iband(find(weight > 0.5)) = ib;
    }
    if (any(iband < 0)) {
      ERROR << "Some g-points not inside a band";
      THROW(1);
    }
    return iband;
  }

  bool is_sw() const {
    return (!solar_irradiance_.empty());
  }

  const Vector& solar_irradiance() const { return solar_irradiance_; }

  void save_g_points(const Vector& wn, const intVector& gp) {
    wavenumber_hr_ = wn;
    g_point_ = gp;
    do_save_g_points_ = true;
  }

  bool read_g_points(Vector& wn, intVector& gp) {
    if (!wavenumber_hr_.empty()) {
      wn = wavenumber_hr_;
      gp = g_point_;
      return true;
    }
    else {
      return false;
    }
  }

  // All the molecules in the CKD model
  std::vector<std::string> molecules;

  // Only the "active" molecules, i.e. those whose coefficients are to
  // be optimized
  std::vector<std::string> active_molecules;

  int ng() { return ng_; }

  int nx() { return x.size(); }

  int ngas() { return single_gas_data_.size(); }

  /// State vector; all the molar absorption coefficient arrays are
  /// linked to parts of this
  Array<1,Real,IsActive> x;

  Vector x_prior;

  bool logarithmic_interpolation = false;
  //bool logarithmic_interpolation = true;
  
  SingleGasData<IsActive>& single_gas(int igas) { return single_gas_data_[igas]; }

  void set_model_id(const std::string mi) { model_id_ = mi; }
  const std::string& model_id() const { return model_id_; }

  /// Calculate Rayleigh molar scattering coefficient in each g point
  void calc_rayleigh_molar_scat(const Vector& ssi_intervals)
  {
    Vector wavenumber_mid = 0.5 * (wavenumber1_ + wavenumber2_);
    // High-res (50 cm-1 resolution) Rayleigh molar scattering
    // coefficient (m2 mol-1)
    Vector rayleigh_molar_scat_hr = rayleigh_molar_scattering_coeff(wavenumber_mid);
    Real ref_surface_pressure = 1.0e5; // Pa
    // Moles of air per m2 of atmospheric column
    Real molar_column = ref_surface_pressure / (ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR);
    Vector optical_depth_hr = molar_column * rayleigh_molar_scat_hr;
    Vector transmission_hr = exp(-optical_depth_hr / REFERENCE_COS_SZA);
    // Weight by solar irradiance
    Vector transmission = (gpoint_fraction_ ** (ssi_intervals*transmission_hr))
                        / (gpoint_fraction_ ** ssi_intervals);
    Vector optical_depth(ng_);
    optical_depth = -log(max(1.0e-14, transmission)) * REFERENCE_COS_SZA;
    rayleigh_molar_scat_ = optical_depth / molar_column;
  }

private:

  /// Is "gas" active, i.e. one for which we want to optimize its
  /// coefficients?
  bool is_active_(const std::string& gas, const std::vector<std::string>& gas_list) {
    bool is_active = false;
    if (IsActive) {
      if (gas_list.empty()) {
	is_active = true;
      }
      else {
	if (std::find(gas_list.begin(), gas_list.end(), gas)
	    != gas_list.end()) {
	  // Gas found
	  is_active = true;
	}
      }
    }
    return is_active;
  }

  /// Molar absorption coefficients and other variables for each gas
  std::vector<SingleGasData<IsActive> > single_gas_data_;

  /// Temperature coordinate variable for Planck function look-up
  /// table (K)
  Vector temperature_planck_;

  /// Planck function at each g point and temperature (W m-2),
  /// dimensioned (temperature,g-point)
  Matrix planck_function_;

  /// Solar irradiance in each g point (W m-2)
  Vector solar_irradiance_;

  /// Rayleigh molar scattering coefficient (m2 mol-1) in each g-point
  Array<1,adept::Real,IsActive> rayleigh_molar_scat_;

  /// Are we retrieving the Rayleigh molar scattering coefficient?
  bool rayleigh_is_active_;

  // Offset to start of state vector where Rayleigh coefficients are
  // stored
  int rayleigh_ix_;

  /// Inverse of diagonal of prior error covariance matrix for
  /// Rayleigh coefficients
  Vector rayleigh_inv_background_;

  /// Natural logarithm of pressure coordinate variable for molar
  /// absorption look-up tables (Pa)
  Vector log_pressure_;

  /// Temperature for molar absorption look-up tables (K), varying
  /// with pressure
  Matrix temperature_;

  /// Bounds of spectral intervals, cm-1, dimensioned (nwav)
  Vector wavenumber1_, wavenumber2_;
  /// Fraction of spectrum contributing to each g-point, dimensioned
  /// (ng,nwav)
  Matrix gpoint_fraction_;

  /// Wavenumber bounds of each band
  Vector wavenumber1_band_, wavenumber2_band_;
  /// Band to which each wavenumber belongs
  intVector band_number_;

  /// Sometimes create_lut changes the number of g points, in which
  /// case it needs to write the full mapping for scale_lut.
  Vector wavenumber_hr_;
  intVector g_point_;
  bool do_save_g_points_ = false;

  /// Number of g points, pressures, temperatures
  int ng_, nt_, np_, nwav_;

  /// Contents of history attribute when reading file
  std::string history_, summary_, config_, model_id_;

};



#endif
