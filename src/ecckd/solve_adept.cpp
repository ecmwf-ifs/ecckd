// solve_adept.cpp - Optimize look-up table using Adept library
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

#include "solve_adept.h"
#include "calc_cost_function_lw.h"
#include "calc_cost_function_sw.h"
#include "Error.h"
#include "Timer.h"

static const Real MIN_X = -1.0e20;

void
calc_total_optical_depth(CkdModel<true>& ckd_model, const LblFluxes& lbl1,
			 aArray3D& optical_depth, bool first_call)
{
  int nprof = lbl1.pressure_hl_.dimension(0);
  int nlev  = lbl1.pressure_hl_.dimension(1)-1;
  optical_depth.resize(nprof,nlev,ckd_model.ng());
  optical_depth = 0.0;

  // Rayleigh scattering
  if (ckd_model.is_sw()) {
    optical_depth += ckd_model.calc_rayleigh_optical_depth(lbl1.pressure_hl_);
  }

  Matrix temperature_fl, p_x_t;
  p_x_t = lbl1.temperature_hl_ * lbl1.pressure_hl_;
  temperature_fl = (p_x_t(__,range(0,end-1)) + p_x_t(__,range(1,end)))
    / (lbl1.pressure_hl_(__,range(0,end-1)) + lbl1.pressure_hl_(__,range(1,end)));


  for (int igas = 0; igas < ckd_model.molecules.size(); ++igas) {
    if (lbl1.gas_mapping(igas) >= 0) {
      if (first_call) {
	LOG << "    Computing CKD optical depth of \"" << ckd_model.molecules[igas] << "\"\n";
      }
      optical_depth += ckd_model.calc_optical_depth(ckd_model.molecules[igas],
						    lbl1.pressure_hl_,
						    temperature_fl,
						    lbl1.vmr_fl_(__,lbl1.gas_mapping(igas),__));
    }
    else {
      if (ckd_model.single_gas(igas).conc_dependence == NONE) {
	if (first_call) {
	  LOG << "    Computing CKD optical depth of background gases in \"" << ckd_model.molecules[igas] << "\"\n";
	}
	optical_depth += ckd_model.calc_optical_depth(ckd_model.molecules[igas],
						      lbl1.pressure_hl_,
						      temperature_fl);
      }
      else {
	if (first_call) {
	  LOG << "    Skipping \"" << ckd_model.molecules[igas] << "\": not present in LBL file\n";
	}
      }
    }
  }
}

Real
calc_cost_function_and_gradient(CkdModel<true>& ckd_model,
				std::vector<LblFluxes>& lbl,
				Vector gradient,
				Real flux_weight, 
				Real flux_profile_weight,
				Real broadband_weight,
				Real spectral_boundary_weight,
				Real erythemal_weight,
				Real negative_od_penalty,
				Real pressure_weight_power,
				Array3D* relative_ckd_flux_dn,
				Array3D* relative_ckd_flux_up)
{
  static bool first_call = true;
  //  data.timer.start(data.rt_id);
  if (first_call) {
    LOG << "  First calculation of cost function and gradient\n";
  }

  ADEPT_ACTIVE_STACK->new_recording();
  aReal cost = 0.0;
  Vector cost_fn_per_band(maxval(lbl[0].iband_per_g)+1);
  cost_fn_per_band = 0.0;

  // Loop over training scenes
  for (int ilbl = 0; ilbl < lbl.size(); ++ilbl) {
    if (first_call) {
      LOG << "  LBL training scene " << ilbl << "\n";
    }
    LblFluxes& lbl1 = lbl[ilbl];
    int nprof = lbl1.pressure_hl_.dimension(0);
    int nlev  = lbl1.pressure_hl_.dimension(1)-1;
    aArray3D optical_depth;
    calc_total_optical_depth(ckd_model, lbl1,
			     optical_depth, first_call);
    int nnegative = count(optical_depth < 0.0);
    if (nnegative > 0) {
      aArray3D penalty(optical_depth.dimensions());
      penalty.where(optical_depth < 0.0) = either_or(optical_depth*optical_depth, 0.0);
      aReal new_cost = negative_od_penalty * sum(penalty);
      cost += new_cost;
      optical_depth.where(optical_depth < 0.0) = 0.0;
      LOG << "  Fixing negative optical depth at " << nnegative << " points: added "
          << new_cost << " to cost function\n";
    }

    // If the pointers relative_ckd_flux_[dn|up] are not NULL, then
    // they point to 3D arrays of fluxes to be subtracted from the CKD
    // calculations.  But since each profile is analyzed in turn, we
    // need to obtain a pointer to the relevant profile, or NULL.
    Matrix* rel_ckd_flux_dn = 0;
    Matrix* rel_ckd_flux_up = 0;
    Matrix  rel_ckd_flux_dn_ref, rel_ckd_flux_up_ref;
    if (relative_ckd_flux_dn) {
      rel_ckd_flux_dn = &rel_ckd_flux_dn_ref;
      rel_ckd_flux_up = &rel_ckd_flux_up_ref;
    }

    // Loop over profiles in one scene
    for (int iprof = 0; iprof < nprof; ++iprof) {
      Vector layer_weight;
      if (pressure_weight_power == 0.5) {
	layer_weight = sqrt(lbl1.pressure_hl_(iprof,range(1,end)))-sqrt(lbl1.pressure_hl_(iprof,range(0,end-1)));
      }
      else if (pressure_weight_power == 1.0) {
	layer_weight = lbl1.pressure_hl_(iprof,range(1,end)) - lbl1.pressure_hl_(iprof,range(0,end-1));
      }
      else {
	layer_weight = pow(lbl1.pressure_hl_(iprof,range(1,end)),pressure_weight_power)
	  -pow(lbl1.pressure_hl_(iprof,range(0,end-1)),pressure_weight_power);
      }
      layer_weight /= sum(layer_weight);

      // Make a soft link to a slice of the relative-to fluxes
      if (relative_ckd_flux_dn) {
	rel_ckd_flux_dn_ref >>= (*relative_ckd_flux_dn)[iprof];
	rel_ckd_flux_up_ref >>= (*relative_ckd_flux_up)[iprof];
      }

      if (!lbl1.is_sw()) {
	Vector spectral_flux_dn_surf, spectral_flux_up_toa;
	if (!lbl1.spectral_flux_dn_surf_.empty()) {
	  spectral_flux_dn_surf >>= lbl1.spectral_flux_dn_surf_(iprof,__);
	  spectral_flux_up_toa  >>= lbl1.spectral_flux_up_toa_(iprof,__);
	}
	cost += calc_cost_function_ckd_lw(lbl1.pressure_hl_(iprof,__),
					  lbl1.planck_hl_(iprof,__,__),
					  lbl1.surf_emissivity_(iprof,__),
					  lbl1.surf_planck_(iprof,__),
					  optical_depth(iprof,__,__),
					  lbl1.spectral_flux_dn_(iprof,__,__),
					  lbl1.spectral_flux_up_(iprof,__,__),
					  lbl1.spectral_heating_rate_(iprof,__,__),
					  spectral_flux_dn_surf,
					  spectral_flux_up_toa,
					  flux_weight, flux_profile_weight, broadband_weight,
					  spectral_boundary_weight,
					  layer_weight, rel_ckd_flux_dn, rel_ckd_flux_up, 
					  lbl1.iband_per_g);
      }
      else {
	//	LOG << "   " << iprof;
	//Real tsi_scaling = sum(lbl1.spectral_flux_dn_(iprof,0,__))
	//  / (lbl1.mu0_(iprof) * sum(ckd_model.solar_irradiance()));
	Real tsi_scaling = lbl1.tsi_ / sum(ckd_model.solar_irradiance());
	Vector spectral_flux_dn_surf, spectral_flux_up_toa;
	if (!lbl1.spectral_flux_dn_surf_.empty()) {
	  spectral_flux_dn_surf >>= lbl1.spectral_flux_dn_surf_(iprof,__);
	  spectral_flux_up_toa  >>= lbl1.spectral_flux_up_toa_(iprof,__);
	}
	Vector spectral_boundary_weights = erythemal_weight * lbl1.erythemal_spectrum_;
	cost += calc_cost_function_ckd_sw(lbl1.mu0_(iprof),
					  lbl1.pressure_hl_(iprof,__),
					  tsi_scaling * ckd_model.solar_irradiance(),
					  lbl1.effective_spectral_albedo_,
					  optical_depth(iprof,__,__),
					  lbl1.spectral_flux_dn_(iprof,__,__),
					  lbl1.spectral_flux_up_(iprof,__,__),
					  lbl1.spectral_heating_rate_(iprof,__,__),
					  spectral_flux_dn_surf,
					  spectral_flux_up_toa,
					  flux_weight, flux_profile_weight, broadband_weight,
					  spectral_boundary_weights,
					  layer_weight, rel_ckd_flux_dn, rel_ckd_flux_up, 
					  lbl1.iband_per_g, cost_fn_per_band);
      }
    } 
  }
  //  data.timer.start(data.autodiff_id);
  cost.set_gradient(1.0);
  ADEPT_ACTIVE_STACK->reverse();
  ckd_model.x.get_gradient(gradient);

  first_call = false;

  //  LOG << cost << " " << maxval(fabs(gradient)) << "\n";
  //  LOG << "  Cost per band = " << cost_fn_per_band << "\n";
  //  data.timer.start(data.minimizer_id);
  return value(cost);
}


struct MyData {
  MyData() {
    minimizer_id = timer.new_activity("minimizer");
    background_id = timer.new_activity("a-priori");
    rt_id = timer.new_activity("radiative transfer");
    //    autodiff_id = timer.new_activity("automatic differentiation");
  }
  CkdModel<true>* ckd_model;
  std::vector<LblFluxes>* lbl;
  Real flux_weight, flux_profile_weight, broadband_weight, prior_error;
  Real spectral_boundary_weight, negative_od_penalty;
  Real pressure_weight_power = 0.5;
  Real erythemal_weight = 0.0;
  Array3D* relative_ckd_flux_dn;
  Array3D* relative_ckd_flux_up;
  Timer timer;
  int minimizer_id, background_id, rt_id, autodiff_id;
};

class CkdOptimizable : public adept::Optimizable {
public:
  virtual Real calc_cost_function(const Vector& x) {
    Vector gradient(x.size());
    return calc_cost_function_gradient(x, gradient);
  }

  virtual Real calc_cost_function_gradient(const Vector& xdata, Vector gradient) {
    aVector& x = data.ckd_model->x;
    for (int ix = 0; ix < x.size(); ix++) {
      if (xdata[ix] > MIN_X) {
	x(ix) = exp(xdata[ix]);
      }
      else {
	x(ix) = 0.0;
      }
    }

    data.timer.start(data.rt_id);
    Real J = calc_cost_function_and_gradient(*(data.ckd_model),
					     *(data.lbl),
					     gradient,
					     data.flux_weight, 
					     data.flux_profile_weight,
					     data.broadband_weight,
					     data.spectral_boundary_weight,
					     data.erythemal_weight,
					     data.negative_od_penalty,
					     data.pressure_weight_power,
					     data.relative_ckd_flux_dn,
					     data.relative_ckd_flux_up);
    data.timer.start(data.background_id);
    // Prior contribution
    //#define OLD_PRIOR 1
#ifdef OLD_PRIOR
    Vector gradient_prior = (1.0/(data.prior_error*data.prior_error)) * (xdata-data.ckd_model->x_prior);
    Real J_prior = 0.5*sum((xdata-data.ckd_model->x_prior) * gradient_prior);

#else
    Vector gradient_prior(data.ckd_model->nx());
    Real J_prior = data.ckd_model->calc_background_cost_function(xdata-data.ckd_model->x_prior, gradient_prior);
#endif
    
    for (int ix = 0; ix < data.ckd_model->nx(); ix++) {
      if (xdata[ix] > MIN_X) {
	gradient[ix] = gradient(ix) * value(x(ix)) + gradient_prior(ix);
      }
      else {
	gradient[ix] = 0.0;
      }
    }
    // Sometimes tiny gradients can cause problems in the minimization
    // - set to zero
    gradient.where(fabs(gradient) < 1.0e-80) = 0.0;
    
    J += J_prior;
    ///  std::cout << J << " " << J_prior << " infinity norm " << maxval(fabs(gradient)) << " " << maxval(fabs(gradient_prior)) << std::endl;

    data.timer.start(data.minimizer_id);
    return J;
  }

  virtual void report_progress(int niter, const Vector& x,
			       Real cost, Real gnorm) {
    LOG << "Iteration " << niter << ": cost function = " << cost 
	<< ", gradient norm = " << gnorm << "\n";
  }

  virtual bool provides_derivative(int order) {
    return (order>=0 && order <= 1);
  }

  MyData data;
};


adept::MinimizerStatus
solve_adept(CkdModel<true>& ckd_model,
	    std::vector<LblFluxes>& lbl,
	    Real flux_weight,
	    Real flux_profile_weight,
	    Real broadband_weight,
	    Real spectral_boundary_weight,
	    Real erythemal_weight,
	    Real prior_error,
	    int max_iterations,
	    Real convergence_criterion,
	    Real negative_od_penalty,
	    Real pressure_weight_power,
	    bool is_bounded,
	    Array3D* relative_ckd_flux_dn,
	    Array3D* relative_ckd_flux_up)
{
  int status=0;

  CkdOptimizable ckd_optimizable;
  adept::Minimizer minimizer(MINIMIZER_ALGORITHM_LIMITED_MEMORY_BFGS);
  minimizer.set_max_iterations(max_iterations);
  minimizer.set_max_step_size(2.0);
  minimizer.set_converged_gradient_norm(convergence_criterion);

  // State vector holds the logarithm of the molar absorption coefficients
  Vector x(ckd_model.nx());

  // Values of zero are held at zero
  x = MIN_X;
  x.where(ckd_model.x > 0.0) = log(ckd_model.x.inactive_link());

  ckd_model.x_prior = x;

  Vector x_min, x_max;
  if (is_bounded) {
    adept::minimizer_initialize_bounds(ckd_model.nx(), x_min, x_max);
    x_min.where(ckd_model.x_min > 0.0) = log(ckd_model.x_min.inactive_link());
    x_max.where(ckd_model.x_max > 0.0) = log(ckd_model.x_max.inactive_link());
    // For g points where the minimum absorption is zero, we set it to
    // a more realistic value, which is twice as far on the negative
    // side (in log space) as x_max is from x, but in case x==x_max,
    // we ensure x_min is no more than x_max-1 (both in log space).
    x_min.where(ckd_model.x_min == 0 && ckd_model.x > 0.0 && ckd_model.x_max > 0.0)
      = min(3.0*x - 2.0*x_max, x_max-1.0);
    intVector bad_loc = find(ckd_model.x_max > 0.0 && x_min>=x_max);
    if (!bad_loc.empty()) {
      WARNING << bad_loc.size() << " bounds on the state variables have x_min>=x_max, starting at index "

	      << bad_loc(0);
      int nbad = bad_loc.size();
      if (nbad > 10) {
	intVector bad_loc2 = bad_loc(range(0,9));
	LOG << "x_min(bad)=" << x_min(bad_loc2) << "...\n";
	LOG << "x_max(bad)=" << x_max(bad_loc2) << "...\n";
	LOG << "ckd_model.x_min(bad)=" << ckd_model.x_min(bad_loc2) << "...\n";
	LOG << "ckd_model.x_prior(bad)=" << ckd_model.x_prior(bad_loc2) << "...\n";
	LOG << "ckd_model.x_max(bad)=" << ckd_model.x_max(bad_loc2) << "...\n";
      }
      else {
	LOG << "x_min(bad)=" << x_min(bad_loc) << "\n";
	LOG << "x_max(bad)=" << x_max(bad_loc) << "\n";
	LOG << "ckd_model.x_min(bad)=" << ckd_model.x_min(bad_loc) << "\n";
	LOG << "ckd_model.x_prior(bad)=" << ckd_model.x_prior(bad_loc) << "\n";
	LOG << "ckd_model.x_max(bad)=" << ckd_model.x_max(bad_loc) << "\n";
      }
      ENDWARNING;
    }
  }

  MyData& data = ckd_optimizable.data;
  data.ckd_model = &ckd_model;
  data.lbl       = &lbl;
  data.flux_weight=flux_weight;
  data.flux_profile_weight=flux_profile_weight;
  data.broadband_weight=broadband_weight;
  data.spectral_boundary_weight = spectral_boundary_weight;
  data.erythemal_weight = erythemal_weight;
  data.prior_error = prior_error;
  data.relative_ckd_flux_dn = relative_ckd_flux_dn;
  data.relative_ckd_flux_up = relative_ckd_flux_up;
  data.negative_od_penalty = negative_od_penalty;
  data.pressure_weight_power = pressure_weight_power;

  LOG << "Optimizing coefficients with Adept LBFGS algorithm: max iterations = "
      << max_iterations << ", convergence criterion = " << convergence_criterion << "\n";
  LOG << "  CKD model interpolation is ";
  if (ckd_model.logarithmic_interpolation) {
    LOG << "LOGARITHMIC\n";
  }
  else {
    LOG << "LINEAR\n";
  }

  data.timer.start(data.minimizer_id);

  // Call Adept minimizer
  if (is_bounded) {
    LOG << "  Minimization is bounded\n";
    LOG << "    number bounded below: " << count(ckd_model.x_min>0)
	<< ", above: " << count(ckd_model.x_max>0)
	<< ", out of non-zero elements: " << count(ckd_model.x > 0) << "\n";
    return minimizer.minimize(ckd_optimizable, x, x_min, x_max);
  }
  else {
    LOG << "  Minimization is unbounded\n";
    return minimizer.minimize(ckd_optimizable, x);
  }
}

