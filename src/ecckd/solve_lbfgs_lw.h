#ifndef SOLVE_LBFGS_LW_H
#define SOLVE_LBFGS_LW H 1

#include <vector>
#include "ckd_model.h"
#include "lbl_fluxes.h"

int solve_lbfgs_lw(CkdModel<true>& ckd_model, std::vector<LblFluxes>& fluxes,
		   Real flux_weight, Real flux_profile_weight, Real broadband_weight, Real prior_error);

std::string lbfgs_status_string(int status);

#endif
