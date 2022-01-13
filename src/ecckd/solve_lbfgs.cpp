// solve_lbfgs.cpp - Optimize look-up table using L-BFGS library
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

#include "lbfgs.h"
#include "solve_lbfgs.h"
#include "solve_adept.h"
#include "calc_cost_function_lw.h"
#include "calc_cost_function_sw.h"
#include "Error.h"
#include "Timer.h"


static const Real MIN_X = -1.0e20;

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
  Array3D* relative_ckd_flux_dn;
  Array3D* relative_ckd_flux_up;
  Timer timer;
  int minimizer_id, background_id, rt_id, autodiff_id;
};

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

  data.timer.start(data.rt_id);
  Vector gradient(data.ckd_model->nx());
  lbfgsfloatval_t J = calc_cost_function_and_gradient(*(data.ckd_model),
						      *(data.lbl),
						      gradient,
						      data.flux_weight, 
						      data.flux_profile_weight,
						      data.broadband_weight,
						      data.spectral_boundary_weight,
						      data.negative_od_penalty,
						      data.relative_ckd_flux_dn,
						      data.relative_ckd_flux_up);
  data.timer.start(data.background_id);
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

  data.timer.start(data.minimizer_id);
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
	    Real spectral_boundary_weight,
	    Real prior_error,
	    int max_iterations,
	    Real convergence_criterion,
	    Real negative_od_penalty,
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
  data.spectral_boundary_weight = spectral_boundary_weight;
  data.prior_error = prior_error;
  data.relative_ckd_flux_dn = relative_ckd_flux_dn;
  data.relative_ckd_flux_up = relative_ckd_flux_up;
  data.negative_od_penalty = negative_od_penalty;

  LOG << "Optimizing coefficients with LBFGS algorithm: max iterations = "
      << max_iterations << ", convergence criterion = " << convergence_criterion << "\n";
  LOG << "  CKD model interpolation is ";
  if (ckd_model.logarithmic_interpolation) {
    LOG << "LOGARITHMIC\n";
  }
  else {
    LOG << "LINEAR\n";
  }

  data.timer.start(data.minimizer_id);
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
