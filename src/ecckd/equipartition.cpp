// equipartition.cpp - Class implementing an algorithm for evenly partitioning a 1D space
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

#include <limits>
#include <valarray>
#include <iostream>

#include "equipartition.h"


// abound <= ascale*abound + bscale*bbound
static void merge_bounds(int ni, ep_real* abound, const ep_real* bbound,
			 ep_real ascale, ep_real bscale) {
  for (int ii = 0; ii <= ni; ++ii) {
    abound[ii] = ascale*abound[ii] + bscale*bbound[ii];
  }
}

static void print_bounds(int ni, ep_real* bounds)
{
  std::cout << "      bounds = [";
  for (int ii = 0; ii <= ni; ++ii) {
    std::cout << " " << bounds[ii];
  }
  std::cout << "]\n";
}

static void print_error(int ni, ep_real* error)
{
  std::cout << "      error  = [";
  for (int ii = 0; ii < ni; ++ii) {
    std::cout << " " << error[ii];
  }
  std::cout << "]\n";
}


/// Return pointer to a string describing the status
const char* ep_status_string(EpStatus status) {
  switch(status) {
  case EP_SUCCESS:
    return "Converged";
  case EP_MAX_ITERATIONS_REACHED:
    return "Maximum iterations reached";
  case EP_RESOLUTION_LIMIT_REACHED:
    return "Resolution limit reached";
  case EP_FAILED_TO_CONVERGE:
    return "Failed to converge";
  case EP_NO_PROGRESS:
    return "No progress made";
  case EP_FAILURE:
    return "Unspecified failure";
  case EP_INPUT_ERROR:
    return "Input error";
  default:
    return "Unknown convergence status";
  }
}

/// Compute statistics from "ni" "error" values: the mean error
/// "mean_error", chi-squared "chi2" (sum of squared differences from
/// the mean), and optionally the fractional standard deviation
/// "frac_std" and the fractional range "frac_range" (only assigned if
/// passed a non-NULL pointer
void ep_stats(int ni, ep_real* error, ep_real& mean_error,
	      ep_real& chi2, ep_real* frac_std, ep_real* frac_range) {
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
    *frac_range = (max_error - min_error) / mean_error;
  }

  chi2 = 0.0;
  for (int ii = 0; ii < ni; ++ii) {
    chi2 += (error[ii]-mean_error)*(error[ii]-mean_error);
  }

  if (frac_std) {
    *frac_std = sqrt(chi2 / ni) / mean_error;
  }
}

/// Print a summary of the result of an equipartition call to standard
/// output
void ep_print_result(EpStatus istatus, int iverbose,
		     int ni, ep_real* bounds, ep_real* error)
{
  ep_real mean_error, chi2;
  ep_real frac_std, frac_range;
  ep_stats(ni, &error[0], mean_error, chi2,
	   &frac_std, &frac_range);

  std::cout << "  Equipartition status: " << ep_status_string(istatus)
	    << "\n      fractional range = " << frac_range
	    << "\n      fractional standard deviation = " << frac_std
	    << std::endl;
  if (iverbose > 1) {
    print_bounds(ni, bounds);
    print_error(ni, error);
  }
}

/// Calculate "cost function" from a list of errors: this is either
/// the "fractional range" (maximum error minus minimum error divided
/// by mean error) or the "fractional standard deviation".
ep_real
Equipartition::cost_function(int ni, ep_real* error) {

  ep_real mean_error = 0.0;
  ep_real min_error = +std::numeric_limits<ep_real>::infinity();
  ep_real max_error = -std::numeric_limits<ep_real>::infinity();
  for (int ii = 0; ii < ni; ++ii) {
    mean_error += error[ii];
    if (error[ii] < min_error) {
      min_error = error[ii];
    }
    if (error[ii] > max_error) {
      max_error = error[ii];
    }
  }
  mean_error /= ni;

  if (do_minimize_frac_range) {
    return (max_error - min_error) / mean_error;
  }
  else {
    ep_real chi2 = 0.0;
    for (int ii = 0; ii < ni; ++ii) {
      chi2 += (error[ii]-mean_error)*(error[ii]-mean_error);
    }
    return sqrt(chi2/ni) / mean_error;
  }
}

EpStatus
Equipartition::line_search(int ni, ep_real* bounds, ep_real* newbounds,
			   ep_real* error) {
  if (!errors_up_to_date) {
    calc_error_all(ni, bounds, error);
    errors_up_to_date = true;
  }
  int iterations_remaining = line_search_max_iterations;
  ep_real start_cost_fn = cost_function(ni, error);
  if (iverbose) {
    std::cout << "      line search" << std::flush; 
  }
  merge_bounds(ni, newbounds, bounds, 0.5, 0.5);
  while (iterations_remaining > 0) {
    calc_error_all(ni, newbounds, error);
    errors_up_to_date = false;
    ep_real cost_fn = cost_function(ni, error);
    if (cost_fn < start_cost_fn) {
      merge_bounds(ni, bounds, newbounds, 0.0, 1.0);
      std::cout << std::endl;
      errors_up_to_date = true;
      return EP_SUCCESS;
    }
    else {
      if (iverbose) {
	std::cout << "." << std::flush;
      }
      merge_bounds(ni, newbounds, bounds, 0.5, 0.5);
    }
    --iterations_remaining;
  }
  if (iverbose) {
    std::cout << "FAILED" << std::endl;
  }
  return EP_NO_PROGRESS;
}

/// Find the optimum bound[1] that partitions two intervals such that
/// error[0] is close to equal error[1]
EpStatus
Equipartition::equipartition_2(ep_real* bounds, ep_real* error)
{
  if (iverbose) {
    std::cout << "      Equipartitioning pair of intervals" << std::flush;
  }
  if (!errors_up_to_date) {
    calc_error_all(2, bounds, error);
    errors_up_to_date = true;
  }

  ep_real bound_left = bounds[0], bound_right = bounds[2];
  ep_real ediff_left, ediff_right;
  ep_real frac_error = 0.5*std::fabs(error[1]-error[0])/(error[0]+error[1]);

  ep_real local_tolerance = partition_tolerance;

  ep_real frac_error_orig = frac_error;
  ep_real newbounds[3] = {bounds[0], bounds[1], bounds[2]};
  ep_real newerror[2]  = {error[0], error[1]};

  int iterations_remaining = partition_max_iterations;

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
      ediff_right = newerror[1]-newerror[0];
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
      ediff_left = newerror[1]-newerror[0];
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
    frac_error = 0.5*std::fabs(ediff)/(newerror[0]+newerror[1]);
    if (frac_error < local_tolerance && frac_error < frac_error_orig) {
      // Converged
      bounds[1] = newbounds[1];
      error[0] = newerror[0];
      error[1] = newerror[1];
      if (iverbose) {
	std::cout << " " << ep_status_string(EP_SUCCESS) << std::endl;
      }
      errors_up_to_date = true;
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

  EpStatus istatus = EP_SUCCESS;

  // Bounds not found; did we reduce the error difference?
  if (frac_error < frac_error_orig) {
    // Yes: save results
    bounds[1] = newbounds[1];
    error[0] = newerror[0];
    error[1] = newerror[1];
    if (iverbose) {
      std::cout << " " << frac_error_orig << "=>" << frac_error << std::flush;
    }
    errors_up_to_date = true;

    if (bound_right-bound_left < resolution) {
      istatus = EP_RESOLUTION_LIMIT_REACHED;
    }
    else if (!iterations_remaining) {
      istatus = EP_MAX_ITERATIONS_REACHED;
    }
  }
  else {
    istatus = EP_NO_PROGRESS;
  }

  if (iverbose) {
    std::cout << " " << ep_status_string(istatus) << std::endl;
  }
  return istatus;
}

/// Partition a 1D space into "ni" intervals bounded by "bounds",
/// which points to a vector of length ni+1.  The user supplies the
/// initial bounds. The end values (0 and ni) are not modified, while
/// the interior bounds (1 to ni-1) are. The error associated with
/// each interval is returned in "error" which should point to a
/// vector of length ni.
EpStatus
Equipartition::equipartition_n(int ni, ep_real* bounds_out, ep_real* error)
{
  if (ni == 2) {
    return equipartition_2(bounds_out, error);
  }

  if (iverbose) {
    std::cout << "  Equipartitioning into " << ni << " intervals, partition tolerance "
	      << partition_tolerance << std::endl;
  }

  EpStatus istatus = EP_SUCCESS;

  int n_shuffle_remaining = partition_max_iterations/2;

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
      std::cout << "    " << iterations_remaining
		<< " iterations remaining" << std::flush;
      if (iverbose > 1) {
	print_bounds(ni, &bounds[0]);
      }
    }

    if (!errors_up_to_date) {
      calc_error_all(ni, &bounds[0], error);
      errors_up_to_date = true;
    }

    if (iverbose > 1) {
      print_error(ni, error);
    }

    ep_real cost_fn = cost_function(ni, error);

    if (iverbose) {
      std::cout << "\n      cost function = "
		<< cost_fn << std::endl;
    }

    if (cost_fn < partition_tolerance) {
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
	if (std::fabs(newbounds[ii]-bounds[ii]) > resolution) {
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
      istatus = EP_FAILED_TO_CONVERGE;
      int nnoprogress = 0;
      if (ni > 2 && n_shuffle_remaining > 0) {
	std::cout << "    Shuffle (" << n_shuffle_remaining << " shuffles remaining)" << std::endl;
	if (n_shuffle_remaining % 2) {
	  // Shuffle low -> high -> low
	  for (int ii = 0; ii < ni-1; ++ii) {
	    EpStatus iistatus = equipartition_2(&bounds[ii], &error[ii]);
	    if (iistatus == EP_NO_PROGRESS) {
	      ++nnoprogress;
	    }
	  }
	  for (int ii = ni-3; ii >= 0; --ii) {
	    EpStatus iistatus = equipartition_2(&bounds[ii], &error[ii]);
	    if (iistatus == EP_NO_PROGRESS) {
	      ++nnoprogress;
	    }
	  }
	}
	else {
	  // Shuffle high -> low -> high
	  for (int ii = ni-2; ii >= 0; --ii) {
	    EpStatus iistatus = equipartition_2(&bounds[ii], &error[ii]);
	    if (iistatus == EP_NO_PROGRESS) {
	      ++nnoprogress;
	    }
	  }
	  for (int ii = 1; ii < ni-1; ++ii) {
	    EpStatus iistatus = equipartition_2(&bounds[ii], &error[ii]);
	    if (iistatus == EP_NO_PROGRESS) {
	      ++nnoprogress;
	    }
	  }
	}
	--n_shuffle_remaining;

	if (cost_function(ni, error) < partition_tolerance) {
	  istatus = EP_SUCCESS;
	  break;
	}
	else if (nnoprogress >= ni*2-3) {
	  istatus = EP_FAILED_TO_CONVERGE;
	}
	else {
	  // Success here means I reduced the fractional range a bit
	  // and should try again with the general minimization
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

  errors_up_to_date = false; // Just in case a subsequent call
			     // provides new bounds...

  return istatus;

}

/// Partition a 1D space bounded by "bound0" and "boundn" into
/// intervals, each with an "error" no more than "target_error". The
/// results are written in terms of the number of intervals needed
/// "ni" and the bounds to those intervals "bounds" (will be resized
/// to length n+1) and "error" (resized to length ni).
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

  if (iverbose) {
    std::cout << "  Working out many intervals are needed for target error of "
	      << target_error << std::endl;
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

  if (iverbose) {
    std::cout << "  " << ni << " intervals needed" << std::endl;
  }

  // Repartition
  errors_up_to_date = true;
  return equipartition_n(error.size(), &bounds[0], &error[0]);
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

  if (iverbose) {
    std::cout << "    Finding next bound below " << bound2;
  }

  if (*error_test_value < 0.0) {
    error_test = calc_error(bound1_test, bound2);
  }
  else {
    error_test = *error_test_value;
  }

  while (iterations_remaining > 0 
	 && (error_test > max_error || error_test < min_error)) {
    if (iverbose) {
      std::cout << "." << std::flush;
    }
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

    if (bound1_low == bound1_high) {
      // No further progress to be made; perhaps we have filled domain
      // with a single interval
      break;
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

  if (iverbose) {
    std::cout << " " << bound1_test << std::endl;
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

  if (iverbose) {
    std::cout << "    Finding next bound above " << bound1;
  }
  if (*error_test_value < 0.0) {
    error_test = calc_error(bound1, bound2_test);
  }
  else {
    error_test = *error_test_value;
  }

  while (iterations_remaining > 0 
	 && (error_test > max_error || error_test < min_error)) {
    if (iverbose) {
      std::cout << "." << std::flush;
    }

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

    if (bound2_low == bound2_high) {
      // No further progress to be made; perhaps we have filled domain
      // with a single interval
      break;
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

  if (iverbose) {
    std::cout << " " << bound2_test << std::endl;
  }

  return bound2_test;
}
