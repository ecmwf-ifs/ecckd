#include "lbfgs.h"
#include "solve_lbfgs_lw.h"
#include "calc_cost_function_lw.h"
#include "Error.h"

static const Real MIN_X = -1.0e20;

struct MyData {
  CkdModel<true>* ckd_model;
  std::vector<LblFluxes>* fluxes;
  Real flux_weight, flux_profile_weight, broadband_weight, prior_error;
};


Real
calc_cost_function_and_gradient(CkdModel<true>& ckd_model,
				std::vector<LblFluxes>& fluxes,
				Vector gradient,
				Real flux_weight, 
				Real flux_profile_weight,
				Real broadband_weight)
{
  ADEPT_ACTIVE_STACK->new_recording();
  aReal cost = 0.0;
  for (int iflux = 0; iflux < fluxes.size(); ++iflux) {
    LblFluxes& data = fluxes[iflux];
    int nprof = data.pressure_hl_.dimension(0);
    int nlev  = data.pressure_hl_.dimension(1)-1;
    aArray3D optical_depth(nprof,nlev,ckd_model.ng());
    optical_depth = 0.0;
    for (int igas = 0; igas < data.ngas(); ++igas) {
      optical_depth += ckd_model.calc_optical_depth(data.molecules_[igas],
						    data.pressure_hl_,
						    data.temperature_hl_,
						    data.vmr_fl_(__,igas,__));
    }
    
    for (int iprof = 0; iprof < nprof; ++iprof) {
      Vector layer_weight = sqrt(data.pressure_hl_(iprof,range(1,end)))-sqrt(data.pressure_hl_(iprof,range(0,end-1)));
      layer_weight /= sum(layer_weight);
      cost += calc_cost_function_ckd_lw(data.pressure_hl_(iprof,__),
					data.planck_hl_(iprof,__,__),
					data.surf_emissivity_(iprof,__),
					data.surf_planck_(iprof,__),
					optical_depth(iprof,__,__),
					data.spectral_flux_dn_(iprof,__,__),
					data.spectral_flux_up_(iprof,__,__),
					data.spectral_heating_rate_(iprof,__,__),
					flux_weight, flux_profile_weight, broadband_weight,
					layer_weight);
    } 
  }
  cost.set_gradient(1.0);
  ADEPT_ACTIVE_STACK->reverse();
  ckd_model.x.get_gradient(gradient);

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
  //LOG << "x=" << x << "\n";

  Vector gradient(data.ckd_model->nx());
  lbfgsfloatval_t J = calc_cost_function_and_gradient(*(data.ckd_model),
						      *(data.fluxes),
						      gradient,
						      data.flux_weight, 
						      data.flux_profile_weight,
						      data.broadband_weight);
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
  //  std::cout << J << " " << J_prior << std::endl; // "\n";

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
solve_lbfgs_lw(CkdModel<true>& ckd_model,
	       std::vector<LblFluxes>& fluxes,
	       Real flux_weight,
	       Real flux_profile_weight,
	       Real broadband_weight,
	       Real prior_error)
{
  int status=0;
  lbfgsfloatval_t fx;
  lbfgs_parameter_t param;
  lbfgs_parameter_init(&param);
  param.epsilon = 0.02;
  param.max_iterations = 1000;
  param.max_step = 2.0;
  param.initial_step_size = 0.5;

  Vector x(ckd_model.nx());
  x = MIN_X;
  x.where(ckd_model.x > 0.0) = log(ckd_model.x.inactive_link());

  ckd_model.x_prior = x;

  MyData data;
  data.ckd_model = &ckd_model;
  data.fluxes    = &fluxes;
  data.flux_weight=flux_weight;
  data.flux_profile_weight=flux_profile_weight;
  data.broadband_weight=broadband_weight;
  data.prior_error = prior_error;

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
