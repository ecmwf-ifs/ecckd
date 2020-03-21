#include <algorithm>
#include <sstream>
#include <cmath>
#include <iostream>

#include "ckd_model.h"
#include "DataFile.h"
#include "OutputDataFile.h"
#include "constants.h"
#include "adept_scalar.h"

using namespace adept;

template<bool IsActive>
void
CkdModel<IsActive>::read(const std::string& file_name, 
			 const std::vector<std::string>& gas_list)
{
  clear();
  LOG << "Reading " << file_name << "\n";
  DataFile file(file_name);
  file.read(temperature_planck_, "temperature_planck");
  file.read(temperature_, "temperature");
  file.read(planck_function_, "planck_function");

  {
    Vector pressure;
    file.read(pressure, "pressure");
    log_pressure_ = log(pressure);
  }

  np_ = log_pressure_.size();
  nt_ = temperature_.dimension(0);
  ng_ = planck_function_.dimension(1);

  std::string molecules_str;
  file.read(molecules_str, DATA_FILE_GLOBAL_SCOPE, "molecules");
  
  int n_gases;
  file.read(n_gases, "n_gases");

  std::string molecule;
  int igas = 0;

  // Count number of state variables
  int nx = 0;

  {
    std::stringstream molecules_s(molecules_str);
    while (std::getline(molecules_s, molecule, ' ')) {
      intVector ndims = file.size(molecule + "_molar_absorption_coeff");
      if (is_active_(std::string(molecule), gas_list)) {
	nx += product(ndims);
      }
    }
  }

  if (IsActive) {
    x.resize(nx);
    LOG << "  State variable size: " << nx << "\n";
  }
  int ix = 0; // Current index to state vector x

  // Read in molar absorption coefficients
  std::stringstream molecules_s(molecules_str);
  while (std::getline(molecules_s, molecule, ' ')) {
    single_gas_data_.push_back(SingleGasData<IsActive>(molecule));
    molecules.push_back(molecule);
    SingleGasData<IsActive>& this_gas = single_gas_data_[igas];

    LOG << "  Reading absorption properties of " << this_gas.Molecule << "\n";
    if (file.exist(molecule + "_vmr")) {
      this_gas.conc_dependence = LUT;
      file.read(this_gas.vmr, molecule + "_vmr");
      int nconc = this_gas.vmr.size();
      if (is_active_(this_gas.molecule, gas_list)) {
	// Link to state vector
	this_gas.is_active = true;
	this_gas.molar_abs_conc.clear();
	this_gas.molar_abs_conc >>= x(range(ix,ix+nconc*nt_*np_*ng_-1)).reshape(dimensions(nconc,nt_,np_,ng_));
	LOG << "    preparing to retrieve array of size " << this_gas.molar_abs_conc.dimensions() << "\n";
	this_gas.ix = ix;
	ix += nconc*nt_*np_*ng_;
      }
      else {
	this_gas.molar_abs_conc.resize(nconc,nt_,np_,ng_);
      }
      for (int iconc = 0; iconc < nconc; ++iconc) {
	Array3 tmp;
	file.read(tmp, molecule + "_molar_absorption_coeff", iconc);
	this_gas.molar_abs_conc(iconc,__,__,__) = tmp;
      }
    }
    else {
      this_gas.conc_dependence = LINEAR;
      Array3 tmp;
      file.read(tmp, molecule + "_molar_absorption_coeff");

      this_gas.molar_abs.clear();
      if (is_active_(this_gas.molecule, gas_list)) {
	this_gas.is_active = true;
	this_gas.molar_abs >>= x(range(ix,ix+nt_*np_*ng_-1)).reshape(dimensions(nt_,np_,ng_));
	LOG << "    preparing to retrieve array of size " << this_gas.molar_abs.dimensions() << "\n";
	this_gas.ix = ix;
	ix += nt_*np_*ng_;
      }
      this_gas.molar_abs = tmp;

      /*
      if (molecule == "co2") {
	LOG << "abs_co2 = " << this_gas.molar_abs(__,__,end) << "\n";
      }
      */

    }
    ++igas;
  }

  if (nx != ix) {
    ERROR << "Mismatch between number of state variables and number of coefficients to optimize";
    THROW(UNEXPECTED_EXCEPTION);
  }

}

template<bool IsActive>
void
CkdModel<IsActive>::write(const std::string& file_name,
			  int argc, const char* argv[], const std::string& config_str)
{
  OutputDataFile file(file_name);

  file.define_dimension("temperature", nt_);
  file.define_dimension("pressure", np_);
  file.define_dimension("g_point", ng_);
  file.define_dimension("temperature_planck", temperature_planck_.size());

  std::string molecule_list;
  for (int igas = 0; igas < ngas(); ++igas) {
    if (igas == 0) {
      molecule_list = single_gas_data_[igas].molecule;
    }
    else {
      molecule_list += " " + single_gas_data_[igas].molecule;
    }
  }

  file.define_variable("n_gases", INT);
  file.write_long_name("Number of gases treated", "n_gases");
  file.write_comment("The gases are listed in the global attribute \"molecules\".", "n_gases");

  file.define_variable("temperature", FLOAT, "temperature", "pressure");
  file.write_long_name("Temperature", "temperature");
  file.write_units("K", "temperature");

  file.define_variable("pressure", FLOAT, "pressure");
  file.write_long_name("Pressure", "pressure");
  file.write_units("Pa", "pressure");

  file.define_variable("temperature_planck", FLOAT, "temperature_planck");
  file.write_long_name("Temperature for Planck function look-up table", "temperature_planck");
  file.write_units("K", "temperature_planck");

  for (int igas = 0; igas < ngas(); ++igas) {
    const SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
    const std::string& Molecule = this_gas.Molecule;
    const std::string& molecule = this_gas.molecule;

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	
      }
      break;
    case LINEAR:
      {
	file.define_variable(molecule + "_" + K_NAME,
			     FLOAT, "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of " + Molecule,
			     molecule + "_" + K_NAME);
	file.write_units("m2 mol-1", molecule + "_" + K_NAME);
      }
      break;
    case LUT:
      {
	file.define_dimension(molecule + "_vmr", this_gas.vmr.size());

	file.define_variable(molecule + "_vmr", FLOAT, molecule + "_vmr");
	file.write_long_name("Volume mixing ratio of " + Molecule + " for look-up table", molecule + "_vmr");
	file.write_units("1", molecule + "_vmr");

	file.define_variable(molecule + "_" + K_NAME,
			     FLOAT, molecule + "_vmr", "temperature", "pressure", "g_point");
	file.write_long_name("Molar absorption coefficient of " + Molecule,
			     molecule + "_" + K_NAME);
	file.write_units("m2 mol-1", molecule + "_" + K_NAME);
      }
    }
  }

  file.define_variable("planck_function", FLOAT, "temperature_planck", "g_point");
  file.write_long_name("Planck function look-up table", "planck_function");
  file.write_units("W m-2", "planck_function");

  std::string title = "Gas optics definition";
  file.write(title, "title");

  file.write(molecule_list, "molecules");
  file.append_history(argc, argv);

  // Write data

  file.end_define_mode();

  file.write(ngas(), "n_gases");
  file.write(eval(exp(log_pressure_)), "pressure");
  file.write(temperature_, "temperature");
  file.write(temperature_planck_, "temperature_planck");

  for (int igas = 0; igas < ngas(); ++igas) {
    SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
    const std::string& molecule = this_gas.molecule;

    switch(this_gas.conc_dependence) {
    case NONE:
      {
	
      }
      break;
    case LINEAR:
      {
	file.write((this_gas.molar_abs).inactive_link(), molecule + "_" + K_NAME);
      }
      break;
    case LUT:
      {
	file.write(this_gas.vmr, molecule + "_vmr");
	for (int iconc = 0; iconc < this_gas.vmr.size(); ++iconc) {
	  file.write((this_gas.molar_abs_conc(iconc,__,__,__)).inactive_link(),
		     molecule + "_" + K_NAME, iconc);
	}
      }
    }
  }

  file.write(planck_function_, "planck_function");

  file.close();

}


/// Create error covariance matrices
template<>
void
CkdModel<true>::create_error_covariances(Real err, Real pressure_corr,
					 Real temperature_corr, Real conc_corr)
{
  for (int igas = 0; igas < ngas(); ++igas) {
    SingleGasData<true>& this_gas = single_gas_data_[igas];
    if (this_gas.conc_dependence == LINEAR) {
      int nx = nt_ * np_;
      LOG << "  Creating " << nx << "x" << nx << " error covariance matrix for " << this_gas.Molecule << "\n";

      intVector t_index_vec(nx), p_index_vec(nx);
      // Soft links:
      intMatrix t_index = t_index_vec.reshape(dimensions(nt_,np_));
      intMatrix p_index = p_index_vec.reshape(dimensions(nt_,np_));

      // SPREAD CONVENTION IS WRONG!!!
      //      t_index = spread<0>(range(0,nt_-1),np_);
      //      p_index = spread<1>(range(0,np_-1),nt_);
      t_index = spread<1>(range(0,nt_-1),np_);
      p_index = spread<0>(range(0,np_-1),nt_);

      SymmMatrix background(nx,nx);

      background = (err*err)*pow(temperature_corr,1.0*abs(spread<0>(t_index_vec,nx)-spread<1>(t_index_vec,nx)));
      background *= pow(pressure_corr,1.0*abs(spread<0>(p_index_vec,nx)-spread<1>(p_index_vec,nx)));

      this_gas.inv_background = inv(background);
    }
    else {
      int nconc = this_gas.vmr.size();
      int nx = nt_ * np_ * nconc;
      LOG << "  Creating " << nx << "x" << nx << " error covariance matrix for " << this_gas.Molecule << "\n";
      intVector t_index_vec(nx), p_index_vec(nx), c_index_vec(nx);
      // Soft links:
      intArray3D t_index = t_index_vec.reshape(dimensions(nconc,nt_,np_));
      intArray3D p_index = p_index_vec.reshape(dimensions(nconc,nt_,np_));
      intArray3D c_index = c_index_vec.reshape(dimensions(nconc,nt_,np_));

      // SPREAD CONVENTION IS WRONG!!!
      //      t_index = spread<2>(spread<0>(range(0,nt_-1),np_),nconc);
      //      p_index = spread<2>(spread<1>(range(0,np_-1),nt_),nconc);
      //      c_index = spread<0>(spread<0>(range(0,nconc-1),nt_),np_);
      t_index = spread<0>(spread<1>(range(0,nt_-1),np_),nconc);
      p_index = spread<0>(spread<0>(range(0,np_-1),nt_),nconc);
      c_index = spread<2>(spread<1>(range(0,nconc-1),nt_),np_);

      Matrix background(nx,nx);

      background = (err*err)*pow(temperature_corr,1.0*abs(spread<0>(t_index_vec,nx)-spread<1>(t_index_vec,nx)));
      background *= pow(pressure_corr,1.0*abs(spread<0>(p_index_vec,nx)-spread<1>(p_index_vec,nx)));
      background *= pow(conc_corr,1.0*abs(spread<0>(c_index_vec,nx)-spread<1>(c_index_vec,nx)));

      this_gas.inv_background = inv(background);
    }
  }

}


/// Return the background contribution to cost function, J, and also
/// the gradient dJ/dx, where delta_x is the difference between the
/// current state and the prior
template<>
Real
CkdModel<true>::calc_background_cost_function(const Vector& delta_x, Vector gradient)
{
  Real cost_fn = 0.0;
  gradient = 0.0;

  for (int igas = 0; igas < ngas(); ++igas) {
    const SingleGasData<true>& this_gas = single_gas_data_[igas];
    if (this_gas.is_active) {
      int ix = this_gas.ix;
      int nx = this_gas.inv_background.dimension(0);
      Vector gradient_local(nx);
      // g-point is fastest varying dimension so need to stride over it
      int nstride = ng_;
      for (int ig = 0; ig < ng_; ++ig) {
	Vector delta_x_local = delta_x(stride(ix,ix+(nx-1)*nstride,nstride));
	//	LOG << "??? " << this_gas.Molecule << " " << delta_x.size() << " " << ix << " " << nx << " " << nstride << " " << ig << " " << delta_x_local.size() << "\n";
	gradient_local = this_gas.inv_background ** delta_x_local;
	cost_fn += 0.5*dot_product(delta_x_local,gradient_local);
	gradient(stride(ix,ix+(nx-1)*nstride,nstride)) += gradient_local;
	++ix;
      }
    }
  }

  return cost_fn;
}



/// Calculate the optical depth at each g-point for multiple
/// atmospheric columns, where the arguments are dimensioned
/// (column,level) and the output dimensioned (column,level,g-point)
template<bool IsActive>
Array<3,Real,IsActive>
CkdModel<IsActive>::calc_optical_depth(int igas,                         ///< Gas number 
				       const Matrix& pressure_hl,        ///< Pressure at half levels (Pa)
				       const Matrix& temperature_fl,     ///< Temperature at full levels (K)
				       const Matrix& vmr_fl) ///< Volume mixing ratio at full levels
{ 
  using std::fmin;
  using std::fmax;
  
  typedef typename scalar<IsActive>::type areal;
  typedef Array<1,Real,IsActive> avector;
  typedef Array<2,Real,IsActive> amatrix;
  typedef Array<3,Real,IsActive> aarray3D;

  // Assume pressure of LUT is evenly spaced in log space
  Real log_p_0 = log_pressure_(0);
  Real d_log_p = log_pressure_(1)-log_pressure_(0); // Spacing

  // Temperature LUT spacing
  Real d_t = temperature_(1,0)-temperature_(0,0);

  SingleGasData<IsActive>& this_gas = single_gas_data_[igas];
  int ncol = pressure_hl.dimension(0);
  int np   = pressure_hl.dimension(1)-1;
  aarray3D od(ncol,np,ng_);
  Real global_weight = 1.0 / (ACCEL_GRAVITY * 0.001 * MOLAR_MASS_DRY_AIR);

  for (int icol = 0; icol < ncol; ++icol) {
    for (int ip = 0; ip < np; ++ip) {
      // Find interpolation points in pressure
      Real log_pressure_fl = log(0.5*(pressure_hl(icol,ip+1)+pressure_hl(icol,ip)));
      Real pindex0 = (log_pressure_fl-log_p_0) / d_log_p;
      pindex0 = fmax(0.0, fmin(pindex0, np_-1.0001));
      int ip0 = static_cast<int>(pindex0);
      Real pweight1 = pindex0 - ip0;
      Real pweight0 = 1.0 - pweight1;
      // Find interpolation points in temperature
      Real t_0 = pweight0*temperature_(0,ip0) + pweight1*temperature_(0,ip0+1);
      Real tindex0 = (temperature_fl(icol,ip)-t_0) / d_t;
      tindex0 = fmax(0.0, fmin(tindex0, nt_-1.0001));
      int it0 = static_cast<int>(tindex0);
      Real tweight1 = tindex0 - it0;
      Real tweight0 = 1.0 - tweight1;

      // Weight
      Real weight = global_weight * vmr_fl(icol,ip) * (pressure_hl(icol,ip+1)-pressure_hl(icol,ip));

      if (this_gas.conc_dependence == LUT) {
	// Find interpolation points in concentration
	Real log_conc = log(vmr_fl(icol,ip));
	Real d_log_c  = log(this_gas.vmr(1)/this_gas.vmr(0));
	Real cindex0  = (log_conc-log(this_gas.vmr(0))) / d_log_c;
	cindex0 = fmax(0.0, fmin(cindex0, this_gas.vmr.size()-1.0001));
	int ic0 = static_cast<int>(cindex0);
	Real cweight1 = cindex0 - ic0;
	Real cweight0 = 1.0 - cweight1;
	// Tri-linear interpolation
	if (!logarithmic_interpolation) {
	  od(icol,ip,__) = weight
	    * (cweight0
	       * (tweight0  * (pweight0 * this_gas.molar_abs_conc(ic0,it0,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0,it0,ip0+1,__))
		  +tweight1 * (pweight0 * this_gas.molar_abs_conc(ic0,it0+1,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0,it0+1,ip0+1,__)))
	       +cweight1
	       * (tweight0  * (pweight0 * this_gas.molar_abs_conc(ic0+1,it0,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0+1,it0,ip0+1,__))
		  +tweight1 * (pweight0 * this_gas.molar_abs_conc(ic0+1,it0+1,ip0,__)
			       +pweight1* this_gas.molar_abs_conc(ic0+1,it0+1,ip0+1,__))));
	}
	else {
	  od(icol,ip,__) = weight
	    * exp(cweight0
		  * (tweight0  * (pweight0 * log(this_gas.molar_abs_conc(ic0,it0,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0,it0,ip0+1,__)))
		     +tweight1 * (pweight0 * log(this_gas.molar_abs_conc(ic0,it0+1,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0,it0+1,ip0+1,__))))
		  +cweight1
		  * (tweight0  * (pweight0 * log(this_gas.molar_abs_conc(ic0+1,it0,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0+1,it0,ip0+1,__)))
		     +tweight1 * (pweight0 * log(this_gas.molar_abs_conc(ic0+1,it0+1,ip0,__))
				  +pweight1* log(this_gas.molar_abs_conc(ic0+1,it0+1,ip0+1,__)))));
	}
      }
      else {
	if (!logarithmic_interpolation) {
	  // Bi-linear interpolation
	  od(icol,ip,__) = weight
	    * (tweight0  * (pweight0 * this_gas.molar_abs(it0,ip0,__)
			    +pweight1* this_gas.molar_abs(it0,ip0+1,__))
	       +tweight1 * (pweight0 * this_gas.molar_abs(it0+1,ip0,__)
			    +pweight1* this_gas.molar_abs(it0+1,ip0+1,__)));
	}
	else {
	  od(icol,ip,__) = weight
	    * exp(tweight0  * (pweight0 * log(this_gas.molar_abs(it0,ip0,__))
			       +pweight1* log(this_gas.molar_abs(it0,ip0+1,__)))
		  +tweight1 * (pweight0 * log(this_gas.molar_abs(it0+1,ip0,__))
			       +pweight1* log(this_gas.molar_abs(it0+1,ip0+1,__))));
	}
      }
    }
  }

  //LOG << "abs " << molecules[igas] << " = " << this_gas.molar_abs(__,__,end) << "\n";

  return od;
}


template<bool IsActive>
Array3
CkdModel<IsActive>::calc_planck_function(const Matrix& temperature_hl)
{
  int ncol = temperature_hl.dimension(0);
  int np   = temperature_hl.dimension(1);
  Array3 planck(ncol,np,ng_);
  for (int icol = 0; icol < ncol; ++icol) {
    planck(icol,__,__) = calc_planck_function(temperature_hl(icol,__));
  }
  return planck;
}


template<bool IsActive>
Matrix
CkdModel<IsActive>::calc_planck_function(const Vector& temperature)
{
  // Planck temperature LUT spacing
  Real d_t = temperature_planck_(1) - temperature_planck_(0);
  Real t0 = temperature_planck_(0);
  int nt = temperature.size();
  Matrix planck(nt,ng_);
  for (int it = 0; it < nt; ++it) {
    Real tindex0 = (temperature(it)-t0) / d_t;
    if (tindex0 >= 0) {
      // Normal interpolation (and extrapolation for high
      // temperatures)
      int it0 = std::min(static_cast<int>(tindex0), temperature_planck_.size()-2);
      Real tweight1 = tindex0 - it0;
      Real tweight0 = 1.0 - tweight1;
      planck(it,__) = tweight0 * planck_function_(it0,__)
   	            + tweight1 * planck_function_(it0+1,__);
    }
    else {
      // Interpolate linearly to zero
      planck(it,__) = (temperature(it)/t0) * planck_function_(0,__);
    }
  }
  return planck;
}

// Instantiate both active and passive versions of CkdModel
template void CkdModel<false>::read(const std::string&,const std::vector<std::string>&);
template void CkdModel<true>::read(const std::string&,const std::vector<std::string>&);
template void CkdModel<false>::write(const std::string&, int argc, const char* argv[], const std::string&);
template void CkdModel<true>::write(const std::string&, int argc, const char* argv[], const std::string&);
template Array<3,Real,false> CkdModel<false>::calc_optical_depth(int,const Matrix&,const Matrix&,const Matrix&);
template Array<3,Real,true> CkdModel<true>::calc_optical_depth(int,const Matrix&,const Matrix&,const Matrix&);
template Array3 CkdModel<false>::calc_planck_function(const Matrix&);
template Array3 CkdModel<true>::calc_planck_function(const Matrix&);
template Matrix CkdModel<false>::calc_planck_function(const Vector&);
template Matrix CkdModel<true>::calc_planck_function(const Vector&);