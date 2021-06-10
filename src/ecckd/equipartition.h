// equipartition.h - Class implementing an algorithm for evenly partitioning a 1D space
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

/// The Equipartition class implements and algorithm for partitioning
/// a 1D space into intervals such that the "error" in each interval
/// is approximately equal.

#ifndef EQUIPARTITION_H
#define EQUIPARTITION_H 1

#include <vector>

#ifdef SINGLE_PRECISION
  typedef float ep_real;
#else
  typedef double ep_real;
#endif

// Possible return values from equipartition functions
typedef enum {
  EP_SUCCESS = 0,
  EP_MAX_ITERATIONS_REACHED,
  EP_FAILED_TO_CONVERGE,
  EP_RESOLUTION_LIMIT_REACHED,
  EP_NO_PROGRESS,
  EP_FAILURE,
  EP_INPUT_ERROR // Error in initial bounds (e.g. not monotonic)
} EpStatus;

/// Return pointer to a string describing the status
const char* ep_status_string(EpStatus status);

/// Compute statistics from "ni" "error" values: the mean error
/// "mean_error", chi-squared "chi2" (sum of squared differences from
/// the mean), and optionally the fractional standard deviation
/// "frac_std" and the fractional range "frac_range" (only assigned if
/// passed a non-NULL pointer
void ep_stats(int ni, ep_real* error, ep_real& mean_error,
	      ep_real& chi2, ep_real* frac_std, ep_real* frac_range);

/// Print a summary of the result of an equipartition call to standard
/// output
void ep_print_result(EpStatus istatus, int iverbose,
		     int ni, ep_real* bounds, ep_real* error);

/// Class for partitioning a 1D space into intervals such that the
/// "error" in each interval is approximately equal.  This class
/// should be inherited by a user class that implements the
/// "calc_error" virtual function. This is suitable for problems where
/// the 1D space is actually discrete and so the error is a
/// non-differentiable function of the bounds, and when the call to
/// "calc_error" is relatively expensive.
class Equipartition {

public:

  /// Partition a 1D space into "ni" intervals bounded by "bounds",
  /// which points to a vector of length ni+1.  The user supplies the
  /// initial bounds. The end values (0 and ni) are not modified,
  /// while the interior bounds (1 to ni-1) are. The error associated
  /// with each interval is returned in "error" which should point to
  /// a vector of length ni.
  EpStatus equipartition_n(int ni, ep_real* bounds, ep_real* error);

  /// Partition a 1D space bounded by "bound0" and "boundn" into
  /// intervals, each with an "error" no more than "target_error". The
  /// results are written in terms of the number of intervals needed
  /// "ni" and the bounds to those intervals "bounds" (will be resized
  /// to length n+1) and "error" (resized to length ni).
  EpStatus equipartition_e(ep_real target_error,
			   ep_real bound0, ep_real boundn,
			   int& ni,
			   std::vector<ep_real>& bounds,
			   std::vector<ep_real>& error);

  /// Calculate "cost function" from a list of errors: this is either
  /// the "fractional range" (maximum error minus minimum error
  /// divided by mean error) or the "fractional standard deviation".
  ep_real cost_function(int ni, ep_real* error);

  /// This function should be implemented by a child class, and should
  /// return the error associated with an interval bounded by "bound1"
  /// and "bound2".
  virtual ep_real calc_error(ep_real bound1, ep_real bound2) = 0;

  /// Calculate the error in all the current intervals
  void calc_error_all(int ni, ep_real* bounds, ep_real* error) {

    if (do_parallel) {
#pragma omp parallel for schedule (dynamic)
      for (int ii = 0; ii < ni; ++ii) {
	error[ii] = calc_error(bounds[ii], bounds[ii+1]);
      }
    }
    else {
      for (int ii = 0; ii < ni; ++ii) {
	error[ii] = calc_error(bounds[ii], bounds[ii+1]);
      }
    }

    if (iverbose) {
      std::cout << "|";
    }

  }

  /// Set verbosity level: 0=quiet, 1=progress, 2=report interim
  /// values for bounds and errors
  void set_verbose(bool v) { iverbose = v ? 1 : 0; }
  void set_verbose(int v) { iverbose = v; }

  /// Set the maximum number of iterations for an equipartition_n call
  void set_partition_max_iterations(int max_it) {
    partition_max_iterations = max_it;
  }

  /// Set the maximum number of iterations for a line search
  void set_line_search_max_iterations(int max_it) {
    line_search_max_iterations = max_it;
  }

  /// Convergence has been achieved when the "cost function" is less
  /// than this number.  The cost function is either the fractional
  /// range of the errors (the maximum error minus the minimum error
  /// all divided by the mean error) or the fractional standard
  /// deviation.
  void set_partition_tolerance(ep_real pt) {
    partition_tolerance = pt;
  }

  /// Use cubic rather than linear interpolation (experience shows
  /// cubic is slower and less likely to converge)
  void set_cubic_interpolation(bool ci) {
    cubic_interpolation = ci;
  }

  /// Specify the resolution of the raw data; differences in the
  /// bounds less than the resolution will not be deemed significant
  void set_resolution(ep_real res) {
    resolution = res;
  }

  /// Do we run "calc_error" in parallel?  Its implementation in child
  /// classes must of course then be thread safe.
  void set_parallel(bool p) {
    do_parallel = p;
  }

  /// Do we minimize the fractional range?  This is the default;
  /// "false" here uses the fractional standard deviation instead.
  void set_minimize_frac_range(bool mfr) {
    do_minimize_frac_range = mfr;
  }

private:
  ep_real next_bound_above(ep_real target_error,
    ep_real bound1, ep_real boundn, 
    ep_real bound2_test, ep_real* error_test = 0);

  ep_real next_bound_below(ep_real target_error,
    ep_real bound0, ep_real bound2,
    ep_real bound1_test, ep_real* error_test = 0);

  /// Find the optimum bound[1] that partitions two intervals such
  /// that error[0] is close to equal error[1]
  EpStatus equipartition_2(ep_real* bounds, ep_real* error);

  EpStatus line_search(int ni, ep_real* bounds, ep_real* newbounds,
		       ep_real* error);

  ep_real next_bound_error_tolerance = 0.05;
  ep_real partition_tolerance = 0.05;
  ep_real resolution = 0.0; // Minimum meaningful distance between bounds
  int next_bound_max_iterations = 20;
  int partition_max_iterations = 20;
  int line_search_max_iterations = 10;

  int iverbose = 0;
  bool cubic_interpolation = false;
  bool do_parallel = false;
  // Do we minimize the fractional standard deviation or the
  // fractional range?
  bool do_minimize_frac_range = true;

  bool errors_up_to_date = false;

};
#endif
