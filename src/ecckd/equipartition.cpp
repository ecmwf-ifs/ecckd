#include <limits>
#include <valarray>
#include "equipartition.h"

EpStatus
Equipartition::equipartition_n(int ni, ep_real* bounds_out, ep_real* error)
{
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
  while (iterations_remaining > 0) {
    calc_error_all(ni, &bounds[0], error);
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
    if ((max_error - min_error) / max_error < partition_tolerance) {
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
      newbounds[inew] = (cum_error[iold+1]-target) * bounds[iold]
	+ (target-cum_error[iold]) * bounds[iold+1];
    }
    bounds = 0.75*newbounds + 0.25*bounds;
    --iterations_remaining;
  }

  for (int ib = 0; ib < ni+1; ++ib) {
    bounds_out[ib] = bounds[ib];
  }
  if (iterations_remaining) {
    return EP_MAX_ITERATIONS_REACHED;
  }
  else {
    return EP_SUCCESS;
  }

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

  // Repartition
  return equipartition_n(error.size(), &bounds[0], &error[0]);

  /*
  // Try two intervals, first split in half
  ep_real bound_mid = 0.5 * (bound0 + boundn);
  ep_real upper_error = calc_error(bound_mid, boundn);
  ep_real lower_error;
  if (upper_error <= target_error) {
    // Upper half has error less than target error; we may need only
    // one interval
    ep_real test_error = calc_error(bound0, boundn);
    if (test_error <= target_error) {
      // One interval is sufficient
      ni = 1;
      bounds.resize(2);
      bounds[0] = bound0;
      bounds[1] = boundn;
      error.resize(1);
      error[0]  = target_error;
      return EP_SUCCESS;
    }
    else {
      // Try lower half
      lower_error = calc_error(bound0, bound_mid);
      if (lower_error <= target_error) {
	// Both halves have an error less than target error, so we
	// only need two intervals
	ni = 2;
	bounds.resize(3);
	bounds[0] = bound0;
	bounds[1] = bound_mid;
	bounds[2] = boundn;
	error.resize(2);
	return equipartition_n(ni, &bounds[0], &error[0]);
      }
      else {
	// Lower half has error above 
      }
    }
  }
  else {
    // Partition from top to middle
    }
  */

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

  if (!error_test_value) {
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

  if (!error_test_value) {
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
