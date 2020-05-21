#include "lbfgs.h"
#include "solve_lbfgs.h"
#include "calc_cost_function_lw.h"
#include "calc_cost_function_sw.h"
#include "Error.h"

static const Real MIN_X = -1.0e20;

struct MyData {
  CkdModel<true>* ckd_model;
  std::vector<LblFluxes>* lbl;
  Real flux_weight, flux_profile_weight, broadband_weight, prior_error;
  Real negative_od_penalty;
  Array3D* relative_ckd_flux_dn;
  Array3D* relative_ckd_flux_up;
};

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
				Real negative_od_penalty,
				Array3D* relative_ckd_flux_dn,
				Array3D* relative_ckd_flux_up)
{
  static bool first_call = true;

  if (first_call) {
    LOG << "  First calculation of cost function and gradient\n";
  }

  ADEPT_ACTIVE_STACK->new_recording();
  aReal cost = 0.0;
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
      Vector layer_weight = sqrt(lbl1.pressure_hl_(iprof,range(1,end)))-sqrt(lbl1.pressure_hl_(iprof,range(0,end-1)));
      layer_weight /= sum(layer_weight);

      // Make a soft link to a slice of the relative-to fluxes
      if (relative_ckd_flux_dn) {
	rel_ckd_flux_dn_ref >>= (*relative_ckd_flux_dn)[iprof];
	rel_ckd_flux_up_ref >>= (*relative_ckd_flux_up)[iprof];
      }

      if (!lbl1.is_sw()) {
	cost += calc_cost_function_ckd_lw(lbl1.pressure_hl_(iprof,__),
					  lbl1.planck_hl_(iprof,__,__),
					  lbl1.surf_emissivity_(iprof,__),
					  lbl1.surf_planck_(iprof,__),
					  optical_depth(iprof,__,__),
					  lbl1.spectral_flux_dn_(iprof,__,__),
					  lbl1.spectral_flux_up_(iprof,__,__),
					  lbl1.spectral_heating_rate_(iprof,__,__),
					  flux_weight, flux_profile_weight, broadband_weight,
					  layer_weight, rel_ckd_flux_dn, rel_ckd_flux_up, 
					  lbl1.iband_per_g);
      }
      else {
	//	LOG << "   " << iprof;
	//Real tsi_scaling = sum(lbl1.spectral_flux_dn_(iprof,0,__))
	//  / (lbl1.mu0_(iprof) * sum(ckd_model.solar_irradiance()));
	Real tsi_scaling = lbl1.tsi_ / sum(ckd_model.solar_irradiance());
	cost += calc_cost_function_ckd_sw(lbl1.mu0_(iprof),
					  lbl1.pressure_hl_(iprof,__),
					  tsi_scaling * ckd_model.solar_irradiance(),
					  optical_depth(iprof,__,__),
					  lbl1.spectral_flux_dn_(iprof,__,__),
					  lbl1.spectral_flux_up_(iprof,__,__),
					  lbl1.spectral_heating_rate_(iprof,__,__),
					  flux_weight, flux_profile_weight, broadband_weight,
					  layer_weight, rel_ckd_flux_dn, rel_ckd_flux_up, 
					  lbl1.iband_per_g);
      }
    } 
  }
  cost.set_gradient(1.0);
  ADEPT_ACTIVE_STACK->reverse();
  ckd_model.x.get_gradient(gradient);

  first_call = false;

  //  LOG << cost << " " << maxval(fabs(gradient)) << "\n";

  return value(cost);
}


extern "C"
lbfgsfloatval_t
calc_cost_function_and_gradient_lbfgs(void *vdata,
				      const lbfgsfloatval_t *xdata,
				      lbfgsfloatval_t *gdata,
				      const int n,
				      const lbfgsfloatval_t step)
{
  MyData& data = *reinterpret_cast<MyData*>(vdata);
  aVector& x = data.ckd_model->x;
  for (int ix = 0; ix < x.size(); ix++) {
    if (xdata[ix] > MIN_X) {
      x(ix) = exp(xdata[ix]);
    }
    else {
      x(ix) = 0.0;
    }
  }

  Vector gradient(data.ckd_model->nx());
  lbfgsfloatval_t J = calc_cost_function_and_gradient(*(data.ckd_model),
						      *(data.lbl),
						      gradient,
						      data.flux_weight, 
						      data.flux_profile_weight,
						      data.broadband_weight,
						      data.negative_od_penalty,
						      data.relative_ckd_flux_dn,
						      data.relative_ckd_flux_up);
  // Prior contribution
  Vector x_data(const_cast<Real *>(xdata), dimensions(x.size()));
  //#define OLD_PRIOR 1
#ifdef OLD_PRIOR
  Vector gradient_prior = (1.0/(data.prior_error*data.prior_error)) * (x_data-data.ckd_model->x_prior);
  Real J_prior = 0.5*sum((x_data-data.ckd_model->x_prior) * gradient_prior);

#else
  Vector gradient_prior(data.ckd_model->nx());
  Real J_prior = data.ckd_model->calc_background_cost_function(x_data-data.ckd_model->x_prior, gradient_prior);
#endif

  for (int ix = 0; ix < data.ckd_model->nx(); ix++) {
    if (xdata[ix] > MIN_X) {
      gdata[ix] = gradient(ix) * value(x(ix)) + gradient_prior(ix);
    }
    else {
      gdata[ix] = 0.0;
    }
  }

  J += J_prior;
  ///  std::cout << J << " " << J_prior << " infinity norm " << maxval(fabs(gradient)) << " " << maxval(fabs(gradient_prior)) << std::endl;

  return J;
}


extern "C"
int
progress_lbfgs(void *data,
	       const lbfgsfloatval_t *x,
	       const lbfgsfloatval_t *g,
	       const lbfgsfloatval_t fx,
	       const lbfgsfloatval_t xnorm,
	       const lbfgsfloatval_t gnorm,
	       const lbfgsfloatval_t step,
	       int n,
	       int k,
	       int ls)
{
  LOG << "Iteration " << k << ": cost function = " << fx << ", gradient norm = " << gnorm << "\n";
  return 0;
}


int
solve_lbfgs(CkdModel<true>& ckd_model,
	    std::vector<LblFluxes>& lbl,
	    Real flux_weight,
	    Real flux_profile_weight,
	    Real broadband_weight,
	    Real prior_error,
	    int max_iterations,
	    Real convergence_criterion,
	    Array3D* relative_ckd_flux_dn,
	    Array3D* relative_ckd_flux_up)
{
  int status=0;
  lbfgsfloatval_t fx;
  lbfgs_parameter_t param;
  lbfgs_parameter_init(&param);
  param.epsilon = convergence_criterion;
  param.max_iterations = max_iterations;
  param.max_step = 2.0;
  param.initial_step_size = 0.5;

  Vector x(ckd_model.nx());
  x = MIN_X;
  x.where(ckd_model.x > 0.0) = log(ckd_model.x.inactive_link());

  ckd_model.x_prior = x;

  MyData data;
  data.ckd_model = &ckd_model;
  data.lbl       = &lbl;
  data.flux_weight=flux_weight;
  data.flux_profile_weight=flux_profile_weight;
  data.broadband_weight=broadband_weight;
  data.prior_error = prior_error;
  data.relative_ckd_flux_dn = relative_ckd_flux_dn;
  data.relative_ckd_flux_up = relative_ckd_flux_up;
  data.negative_od_penalty = 1.0e5;

  LOG << "Optimizing coefficients with LBFGS algorithm: max iterations = "
      << max_iterations << ", convergence criterion = " << convergence_criterion << "\n";
  LOG << "  CKD model interpolation is ";
  if (ckd_model.logarithmic_interpolation) {
    LOG << "LOGARITHMIC\n";
  }
  else {
    LOG << "LINEAR\n";
  }

  status = lbfgs(ckd_model.nx(), x.data(), &fx,
		 calc_cost_function_and_gradient_lbfgs,
		 progress_lbfgs, &data, &param);

  return status;
}

std::string
lbfgs_status_string(int status)
{
  if (status == LBFGS_SUCCESS) {
    return "Converged";
  }
  else if (status == LBFGS_ALREADY_MINIMIZED) {
    return "Already minimized";
  }
  else if (status < LBFGSERR_OUTOFINTERVAL) {
    return "Out of interval";
  }
  else if (status == LBFGSERR_MAXIMUMLINESEARCH) {
    return "Maximum line search";
  }
  else if (status == LBFGSERR_MAXIMUMITERATION) {
    return "Maximum iterations reached";
  }
  else {
    return "Failed to converge";
  }
}
