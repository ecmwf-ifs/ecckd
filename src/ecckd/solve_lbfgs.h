// solve_lbfgs.h - Optimize look-up table using L-BFGS library
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

#ifndef SOLVE_LBFGS_H
#define SOLVE_LBFGS H 1

#include <vector>
#include "ckd_model.h"
#include "lbl_fluxes.h"
void calc_total_optical_depth(CkdModel<true>& ckd_model, const LblFluxes& lbl1,
			      aArray3D& optical_depth, bool first_call = false);

int solve_lbfgs(CkdModel<true>& ckd_model, std::vector<LblFluxes>& lbl,
		Real flux_weight, Real flux_profile_weight, Real broadband_weight,
		Real spectral_boundary_weight, Real prior_error,
		int max_iterations, Real convergence_criterion,
		Array3* relative_ckd_flux_dn = 0, Array3* relative_ckd_flux_up = 0);

std::string lbfgs_status_string(int status);

#endif
