#include <limits>
#include <valarray>
#include <iostream>

#include "equipartition.h"

ep_real ep_cost_function(int ni, ep_real* error) {

  ep_real mean_error = 0.0;
  for (int ii = 0; ii < ni; ++ii) {
    mean_error += error[ii];
  }
  mean_error /= ni;

  ep_real cost_fn = 0.0;
  for (int ii = 0; ii < ni; ++ii) {
    cost_fn += (error[ii]-mean_error)*(error[ii]-mean_error);
  }

  return cost_fn;
}

void ep_stats(int ni, ep_real* error, ep_real& mean_error,
	      ep_real& cost_fn, ep_real* frac_std, ep_real* frac_range) {
  mean_error = 0.0;
  for (int ii = 0; ii < ni; ++ii) {
    mean_error += error[ii];
  }
  mean_error /= ni;

  ep_real min_error = +std::numeric_limits<ep_real>::infinity();
  ep_real max_error = -std::numeric_limits<ep_real>::infinity();
  for (int ii = 0; ii < ni; ++ii) {
    if (error[ii] < min_error) {
      min_error = error[ii];
    }
    if (error[ii] > max_error) {
      max_error = error[ii];
    }
  }

  if (frac_range) {
    *frac_range = 0.5 * (max_error - min_error) / mean_error;
  }

  cost_fn = 0.0;
  for (int ii = 0; ii < ni; ++ii) {
    cost_fn += (error[ii]-mean_error)*(error[ii]-mean_error);
  }

  if (frac_std) {
    *frac_std = sqrt(cost_fn / ni) / mean_error;
  }
  
}

static
int loc_min_error(int ni, const ep_real* error) {
  ep_real min_error = +std::numeric_limits<ep_real>::infinity();
  int imin = -1;
  for (int ii = 0; ii < ni; ++ii) {
    if (error[ii] < min_error) {
      min_error = error[ii];
      imin = ii;
    }
  }
  return imin;
}

static
int loc_max_error(int ni, const ep_real* error) {
  ep_real max_error = -std::numeric_limits<ep_real>::infinity();
  int imax = -1;
  for (int ii = 0; ii < ni; ++ii) {
    if (error[ii] > max_error) {
      max_error = error[ii];
      imax = ii;
    }
  }
  return imax;
}

// abound <= ascale*abound + bscale*bbound
static
void merge_bounds(int ni, ep_real* abound, ep_real* bbound,
		  ep_real ascale, ep_real bscale) {
  for (int ii = 0; ii <= ni; ++ii) {
    abound[ii] = ascale*abound[ii] + bscale*bbound[ii];
  }
}


EpStatus
Equipartition::line_search(int ni, ep_real* bounds, ep_real* newbounds,
			   ep_real* error) {
  int iterations_remaining = line_search_max_iterations;
  ep_real start_cost_fn = ep_cost_function(ni, error);
  if (iverbose) {
    std::cout << "Line search: " << start_cost_fn;
  }
  merge_bounds(ni, newbounds, bounds, 0.5, 0.5);
  while (iterations_remaining > 0) {
    calc_error_all(ni, newbounds, error);
    ep_real cost_fn = ep_cost_function(ni, error);
    if (cost_fn < start_cost_fn) {
      merge_bounds(ni, bounds, newbounds, 0.0, 1.0);
      if (iverbose) {
	std::cout << " " << cost_fn << "\n";
      }
      return EP_SUCCESS;
    }
    else {
      if (iverbose) {
	std::cout << ".";
      }
      merge_bounds(ni, newbounds, bounds, 0.5, 0.5);
    }
    --iterations_remaining;
  }
  if (iverbose) {
    std::cout << " FAILED\n";
  }
  return EP_FAILED_TO_CONVERGE;
}

EpStatus
Equipartition::equipartition_2(ep_real* bounds, ep_real* error)
{
  ep_real bound_left = bounds[0], bound_right = bounds[2];
  ep_real ediff_left, ediff_right;
  ep_real frac_error = 0.5*abs(error[1]-error[0])/(error[0]+error[1]);

  ep_real local_tolerance = partition_tolerance;

  ep_real frac_error_orig = frac_error;
  ep_real newbounds[3] = {bounds[0], bounds[1], bounds[2]};
  ep_real newerror[2]  = {error[0], error[1]};

  int iterations_remaining = 10;
  if (error[0] > error[1]) {
    // Found an upper (right) limit on the middle bound; now find a
    // lower (left) limit
    bound_right = bounds[1];
    ediff_right = error[1]-error[0];  
    while (iterations_remaining) { // && frac_error > local_tolerance) {
      newbounds[1] = (-ediff_right*newbounds[0] + (newerror[0]+ediff_right)*newbounds[1])
	/ newerror[0];
      calc_error_all(2, newbounds, newerror);
      if (newerror[0] < newerror[1]) {
	// Found left limit
	bound_left = newbounds[1];
	ediff_left = newerror[1]-newerror[0];
	break;
      }
      --iterations_remaining;
    }
  }
  else {
    // Found a lower (left) limit on the middle bound; now find an
    // upper (right) limit
    bound_left = bounds[1];
    ediff_left = error[1]-error[0];
    while (iterations_remaining) { // && frac_error > local_tolerance) {
      newbounds[1] = (ediff_left*newbounds[2] + (newerror[1]-ediff_left)*newbounds[1])
	/ newerror[1];
      calc_error_all(2, newbounds, newerror);
      if (newerror[0] > newerror[1]) {
	// Found right limit
	bound_right = newbounds[1];
	ediff_right = newerror[1]-newerror[0];
	break;
      }
      --iterations_remaining;
    }
  }

  bool no_progress = false;
  ep_real prev_frac_error = frac_error;

  // We have bounds on both sides
  while (iterations_remaining) {
    // Remember ediff_right is negative
    if (no_progress) {
      newbounds[1] = 0.5 * (bound_right + bound_left);
    }
    else {
      newbounds[1] = (ediff_left*bound_right - ediff_right*bound_left) / (ediff_left-ediff_right);
    }
    calc_error_all(2, newbounds, newerror);
    ep_real ediff = newerror[1]-newerror[0];
    frac_error = 0.5*abs(ediff)/(newerror[0]+newerror[1]);
    if (frac_error < local_tolerance && frac_error < frac_error_orig) {
      // Converged
      bounds[1] = newbounds[1];
      error[0] = newerror[0];
      error[1] = newerror[1];
      return EP_SUCCESS;
    }
    else if (frac_error == prev_frac_error) {
      if (no_progress) {
	break;
      }
      else {
	no_progress = true;
      }
    }
    if (ediff < 0) {
      ediff_right = ediff;
      bound_right = newbounds[1];
    }
    else {
      ediff_left = ediff;
      bound_left = newbounds[1];
    }
    prev_frac_error = frac_error;
    --iterations_remaining;
  }

  // Bounds not found; did we reduce the error difference?
  if (frac_error < frac_error_orig) {
    // Yes: save results
    bounds[1] = newbounds[1];
    error[0] = newerror[0];
    error[1] = newerror[1];
  }
  if (bound_right-bound_left < resolution) {
    return EP_RESOLUTION_LIMIT_REACHED;
  }
  else {
    return EP_FAILED_TO_CONVERGE;
  }

}

EpStatus
Equipartition::equipartition_n(int ni, ep_real* bounds_out, ep_real* error)
{
  EpStatus istatus = EP_SUCCESS;

  int n_shuffle_remaining = 5;

  // Check initial bounds are monotonic and increasing
  for (int ib = 0; ib < ni; ++ib) {
    if (bounds_out[ib+1] <= bounds_out[ib]) {
      return EP_INPUT_ERROR;
    }
  }
  std::valarray<ep_real> bounds(ni+1);

  for (int ib = 0; ib < ni+1; ++ib) {
    bounds[ib] = bounds_out[ib];
  }

  int iterations_remaining = partition_max_iterations;

  if (ni == 2) {
    iterations_remaining = 5;
  }

  while (iterations_remaining > 0) {
    if (iverbose) {
      std::cout << "Iterations remaining = " << iterations_remaining << "\n";
      print_bounds(ni, &bounds[0]);
    }
    calc_error_all(ni, &bounds[0], error);
    if (iverbose) {
      std::cout << "  error = ";
      print_error(ni, error);
      std::cout << "\n  cost = " << ep_cost_function(ni, error) << "\n";
    }
    else if (dump_error_iteration) {
      print_error(ni, error);
      std::cout << "\n";
    }

    ep_real mean_error, cost_fn;
    ep_real frac_std, frac_range;

    ep_stats(ni, error, mean_error, cost_fn,
	     &frac_std, &frac_range);
    if (frac_range < partition_tolerance) {
      break;
    }
    std::valarray<ep_real> cum_error(ni+1);
    cum_error[0] = 0.0;
    for (int ii = 0; ii < ni; ++ii) {
      cum_error[ii+1] = cum_error[ii] + error[ii];
    }
    ep_real target_error = cum_error[ni] / ni;

    std::valarray<ep_real> newbounds(ni+1);
    newbounds[0]  = bounds[0];
    newbounds[ni] = bounds[ni];
    int iold = 0;
    for (int inew = 1; inew < ni; ++inew) {
      ep_real target = target_error * inew;
      while (cum_error[iold+1] < target) {
	++iold;
      }
      if (cubic_interpolation) {
	ep_real u = (target - cum_error[iold]) / (cum_error[iold+1]-cum_error[iold]);
	ep_real u2 = u*u;
	ep_real u3 = u2*u;
	// Gradient across the interval
	ep_real grad = (bounds[iold+1]-bounds[iold]) / (cum_error[iold+1]-cum_error[iold]);
	ep_real grad0; // Gradient at start of interval
	if (iold == 0) {
	  // We are in first interval
	  grad0 = grad;
	}
	else {
	  ep_real prior_grad = (bounds[iold]-bounds[iold-1])
	    / (cum_error[iold]-cum_error[iold-1]);
	  grad0 = 0.5*(prior_grad+grad);
	  grad0 = (bounds[iold+1]-bounds[iold-1])
	    / (cum_error[iold+1]-cum_error[iold-1]);
	}
	ep_real grad1; // Gradient at end of interval
	if (iold == ni-1) {
	  // We are in last interval
	  grad1 = grad;
	}
	else {
	  ep_real next_grad = (bounds[iold+2]-bounds[iold+1])
	    / (cum_error[iold+2]-cum_error[iold+1]);
	  grad1 = 0.5*(grad+next_grad);
	  grad1 = (bounds[iold+2]-bounds[iold])
	    / (cum_error[iold+2]-cum_error[iold]);
	}
	/*
	// Check for monotonicity
	ep_real alpha = grad0 / grad;
	ep_real beta = grad1 / grad;
	std::cout << "CHECK = " << alpha*alpha + beta*beta;
	if (alpha*alpha + beta*beta > 9) {
	  std::cout << " ***";
	}
	std::cout << "\n";
	*/
	newbounds[inew] = (2.0*u3-3.0*u2+1) * bounds[iold]
	  + (u3-2.0*u2+u) * grad0
	  + (-2.0*u3+3.0*u2) * bounds[iold+1]
	  + (u3-u2) * grad1;
      }
      else { // Linear interpolation
	newbounds[inew] = ((cum_error[iold+1]-target) * bounds[iold]
			   + (target-cum_error[iold]) * bounds[iold+1])
	  / (cum_error[iold+1]-cum_error[iold]);
      }
    }
    // Check if resolution limit reached: all newbounds are within
    // "resolution" of old bounds
    if (resolution > 0.0) {
      bool found_significant_jump = false;
      for (int ii = 1; ii < ni; ++ii) {
	if (std::abs(newbounds[ii]-bounds[ii]) > resolution) {
	  found_significant_jump = true;
	  break;
	}
      }
      if (!found_significant_jump) {
	for (int ib = 0; ib < ni+1; ++ib) {
	  bounds_out[ib] = bounds[ib];
	}
	return EP_RESOLUTION_LIMIT_REACHED;
      }
    }

    EpStatus ls_status = line_search(ni, &bounds[0], &newbounds[0], error);
    if (ls_status != EP_SUCCESS) {
      istatus = ls_status;
      int nresreached = 0;
      if (ni > 3 && n_shuffle_remaining > 0) {
	for (int ii = 0; ii < ni-1; ++ii) {
	  EpStatus iistatus = equipartition_2(&bounds[ii], &error[ii]);
	  if (iistatus = EP_RESOLUTION_LIMIT_REACHED) {
	    ++nresreached;
	  }
	}
	for (int ii = ni-3; ii >= 0; --ii) {
	  EpStatus iistatus = equipartition_2(&bounds[ii], &error[ii]);
	  if (iistatus = EP_RESOLUTION_LIMIT_REACHED) {
	    ++nresreached;
	  }
	}
	--n_shuffle_remaining;

	ep_stats(ni, error, mean_error, cost_fn,
		 &frac_std, &frac_range);
	if (frac_range < partition_tolerance) {
	  istatus = EP_SUCCESS;
	  break;
	}
	else if (nresreached >= ni*2-3) {
	  istatus = EP_RESOLUTION_LIMIT_REACHED;
	}
	else {
	  istatus = EP_SUCCESS;
	}
      }
      if (istatus != EP_SUCCESS) {
      	break;
      }

    }
    --iterations_remaining;
  }

  for (int ib = 0; ib < ni+1; ++ib) {
    bounds_out[ib] = bounds[ib];
  }
  if (iterations_remaining == 0) {
    istatus = EP_MAX_ITERATIONS_REACHED;
  }

  return istatus;

}

EpStatus
Equipartition::equipartition_e(ep_real target_error,
			       ep_real bound0, ep_real boundn,
			       int& ni,
			       std::vector<ep_real>& bounds,
			       std::vector<ep_real>& error)
{
  // Check bounds increasing
  if (boundn <= bound0) {
    return EP_INPUT_ERROR;
  }

  // Find uppermost interval
  ep_real upper_error = -1.0;
  ep_real upper_bound = next_bound_below(target_error,
					 bound0, boundn,
					 0.05*bound0+0.95*boundn,
					 &upper_error);
  if (upper_bound == bound0) {
    // One interval is sufficient
    ni = 1;
    bounds.resize(2);
    bounds[0] = bound0;
    bounds[1] = boundn;
    error.resize(1);
    error[0]  = upper_error;
    return EP_SUCCESS;
  }

  bounds.resize(1);
  error.clear();
  bounds[0] = bound0;

  // Add lower intervals up to upper interval
  int iint = 0;
  while (bounds[iint] < upper_bound) {
    error.push_back(-1.0);
    bounds.push_back(next_bound_above(target_error,
				      bounds[iint], upper_bound,
				      0.25*bounds[iint]+0.75*upper_bound,
				      &error[iint]));
    ++iint;
  }
  // Push the upper interval onto the end
  error.push_back(upper_error);
  bounds.push_back(boundn);

  ni = error.size();

  // Repartition
  return equipartition_n(error.size(), &bounds[0], &error[0]);

}

void
Equipartition::print_bounds(int ni, ep_real* bounds)
{
  std::cout << "  bounds = [";
  for (int ii = 0; ii <= ni; ++ii) {
    std::cout << " " << bounds[ii];
  }
  std::cout << "]\n";
}

void
Equipartition::print_error(int ni, ep_real* error)
{
  for (int ii = 0; ii < ni; ++ii) {
    std::cout << " " << error[ii];
  }
}

ep_real
Equipartition::next_bound_below(ep_real target_error,
				ep_real bound0, ep_real bound2,
				ep_real bound1_test, ep_real* error_test_value)
{
  ep_real max_error = target_error;
  ep_real min_error = target_error * (1.0 - next_bound_error_tolerance);
  ep_real bound1_low  = bound0;
  ep_real bound1_high = bound2;
  ep_real error_low  = -1.0;
  ep_real error_high = 0.0;
  ep_real error_test;
  int iterations_remaining = next_bound_max_iterations;

  if (*error_test_value < 0.0) {
    error_test = calc_error(bound1_test, bound2);
  }
  else {
    error_test = *error_test_value;
  }

  while (iterations_remaining > 0 
	 && (error_test > max_error || error_test < min_error)) {
    if (error_test > target_error) {
      // Found a lower possible value for bound1
      bound1_low = bound1_test;
      error_low  = error_test;
    }
    else {
      // Found an upper possible value for bound1
      bound1_high = bound1_test;
      error_high  = error_test;
    }
    
    if (error_low > 0.0) {
      // Bracketed: interpolate
      bound1_test = ((target_error-error_high) * bound1_low
		     + (error_low-target_error)* bound1_high)
	/ (error_low - error_high);
      if (error_high == 0.0) {
	// Not yet found a true upper possible value for bound1 so
	// relationship unlikely to be linear: move closer to bound2
	bound1_test = 0.5 * (bound1_test + bound1_high);
      }
      else if (error_test < min_error
	       && error_low > 2.0 * max_error) {
	// Lower possible bound1 too far away to be a likely useful
	// interpolation point
	bound1_test = 0.75*bound1_test + 0.25*bound1_low;
      }
    }
    else {
      // Not bracketed: extrapolate, but not too far
      bound1_test = std::max(bound1_low, bound1_high
	     -0.5*target_error*(bound2-bound1_high) / error_high);
    }

    error_test = calc_error(bound1_test, bound2);
    --iterations_remaining;
  }

  if (error_test_value) {
    *error_test_value = error_test;
  }

  return bound1_test;
}

ep_real
Equipartition::next_bound_above(ep_real target_error,
				ep_real bound1, ep_real boundn,
				ep_real bound2_test, ep_real* error_test_value)
{
  ep_real max_error = target_error;
  ep_real min_error = target_error * (1.0 - next_bound_error_tolerance);
  ep_real bound2_low  = bound1;
  ep_real bound2_high = boundn;
  ep_real error_low  = 0.0;
  ep_real error_high = -1.0;
  ep_real error_test;
  int iterations_remaining = next_bound_max_iterations;

  if (*error_test_value < 0.0) {
    error_test = calc_error(bound1, bound2_test);
  }
  else {
    error_test = *error_test_value;
  }

  while (iterations_remaining > 0 
	 && (error_test > max_error || error_test < min_error)) {
    if (error_test > target_error) {
      // Found an upper possible value for bound2
      bound2_high = bound2_test;
      error_high  = error_test;
    }
    else {
      // Found a lower possible value for bound2
      bound2_low = bound2_test;
      error_low  = error_test;
    }
    
    if (error_high > 0.0) {
      // Bracketed: interpolate
      bound2_test = ((target_error-error_low) * bound2_high
		     + (error_high-target_error)* bound2_low)
	/ (error_high - error_low);
      if (error_low == 0.0) {
	// Not yet found a true lower possible value for bound2 so
	// relationship unlikely to be linear: move closer to bound1
	bound2_test = 0.5 * (bound2_test + bound2_low);
      }
      else if (error_test < min_error
	       && error_low > 2.0 * max_error) {
	// Upper possible bound2 too far away to be a likely useful
	// interpolation point
	bound2_test = 0.75*bound2_test + 0.25*bound2_high;
      }
    }
    else {
      // Not bracketed: extrapolate, but not too far
      bound2_test = std::max(bound2_high, bound2_high
		     -0.5*target_error*(bound2_low-bound1) / error_low);
    }

    error_test = calc_error(bound1, bound2_test);
    --iterations_remaining;
  }

  if (error_test_value) {
    *error_test_value = error_test;
  }

  return bound2_test;
}
