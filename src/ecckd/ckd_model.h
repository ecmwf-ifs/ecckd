#ifndef CKD_MODEL_H
#define CKD_MODEL_H 1

#include <vector>
#include <string>

#include <adept_arrays.h>

#include "Error.h"

using namespace adept;

// 0=none, 1=linear, 2=look-up table
typedef enum {
  NONE = 0,
  LINEAR,
  LUT
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

  // Volume mixing ratio coordinate variable
  Vector vmr;

  // Inverse of the background error covariance matrix for one g point
  SymmMatrix inv_background;

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

  CkdModel(const std::string& file_name,
	   const std::vector<std::string>& gas_list = std::vector<std::string>()) {
    read(file_name, gas_list);
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
			    const Matrix& temperature_hl,     ///< Temperature at half levels (K)
			    const Matrix& vmr_fl = Matrix()); ///< Volume mixing ratio at full levels

  /// As above but provide the name of the gas
  Array<3,Real,IsActive> calc_optical_depth(const std::string& gas, ///< Gas name, lower case
			    const Matrix& pressure_hl,        ///< Pressure at half levels (Pa)
			    const Matrix& temperature_hl,     ///< Temperature at half levels (K)
			    const Matrix& vmr_fl = Matrix()) ///< Volume mixing ratio at full levels
    {
      auto itgas = std::find(molecules.begin(), molecules.end(), gas);
      if (itgas == molecules.end()) {
	ERROR << "CKD model does not contain " << gas;
	THROW(UNEXPECTED_EXCEPTION);
      }
      else {
	int igas = std::distance(molecules.begin(), itgas);
	return calc_optical_depth(igas, pressure_hl, temperature_hl, vmr_fl);
      }
    }

  /// Create error covariance matrices
  void create_error_covariances(Real err, Real pressure_corr, Real temperature_corr, Real conc_corr);

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

  std::vector<std::string> molecules;

  int ng() { return ng_; }

  int nx() { return x.size(); }

  int ngas() { return single_gas_data_.size(); }

  /// State vector; all the molar absorption coefficient arrays are
  /// linked to parts of this
  Array<1,Real,IsActive> x;

  Vector x_prior;

  bool logarithmic_interpolation = false;

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

  /// Temperature coordinate variable for Planck function look-up
  /// table (K)
  Vector temperature_planck_;

  /// Planck function at each g point and temperature (W m-2),
  /// dimensioned (temperature,g-point)
  Matrix planck_function_;

  /// Natural logarithm of pressure coordinate variable for molar
  /// absorption look-up tables (Pa)
  Vector log_pressure_;

  /// Temperature for molar absorption look-up tables (K), varying
  /// with pressure
  Matrix temperature_;

  /// Molar absorption coefficients and other variables for each gas
  std::vector<SingleGasData<IsActive> > single_gas_data_;

  /// Number of g points, pressures, temperatures
  int ng_, nt_, np_;

};



#endif
