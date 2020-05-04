#ifndef SOLVE_LBFGS_H
#define SOLVE_LBFGS H 1

#include <vector>
#include "ckd_model.h"
#include "lbl_fluxes.h"
void calc_total_optical_depth(CkdModel<true>& ckd_model, const LblFluxes& lbl1,
			      aArray3D& optical_depth, bool first_call = false);

int solve_lbfgs(CkdModel<true>& ckd_model, std::vector<LblFluxes>& lbl,
		Real flux_weight, Real flux_profile_weight, Real broadband_weight, Real prior_error,
		Real convergence_criterion,
		Array3* relative_ckd_flux_dn = 0, Array3* relative_ckd_flux_up = 0);

std::string lbfgs_status_string(int status);

#endif
