#ifndef EQUIPARTITION_H
#define EQUIPARTITION_H 1

#include <vector>

#ifdef SINGLE_PRECISION
  typedef float ep_real;
#else
  typedef double ep_real;
#endif

typedef enum {
  EP_SUCCESS = 0,
  EP_MAX_ITERATIONS_REACHED,
  EP_FAILED_TO_CONVERGE,
  EP_RESOLUTION_LIMIT_REACHED,
  EP_FAILURE,
  EP_INPUT_ERROR // Error in initial bounds (e.g. not monotonic)
} EpStatus;

inline
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
  case EP_FAILURE:
    return "Failure";
  case EP_INPUT_ERROR:
    return "Input error";
  default:
    return "Unknown convergence status";
  }
}

ep_real ep_cost_function(int ni, ep_real* error);
void ep_stats(int ni, ep_real* error, ep_real& mean_error,
	      ep_real& cost_fn, ep_real* frac_std, ep_real* frac_range);

/// Class for partitioning a 1D space such that the "error" in each
/// space is approximately equal.  This class should be inherited by a
/// user class that implements the "calc_error" virtual function. This
/// is suitable for problems where the 1D space is actually discrete
/// and so the error is a non-differentiable function of the bounds,
/// and when the call to "calc_error" is relatively expensive.
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

  /// Find the optimum bound[1] that partitions two intervals such
  /// that error[0] is close to equal error[1]
  EpStatus equipartition_2(ep_real* bounds, ep_real* error);

  /// This function should be implemented by a child class, and should
  /// return the error associated with an interval bounded by "bound1"
  /// and "bound2".
  virtual ep_real calc_error(ep_real bound1, ep_real bound2) = 0;

  /// Calculate the error in all the current intervals
  void calc_error_all(int ni, ep_real* bounds, ep_real* error) {
    for (int ii = 0; ii < ni; ++ii) {
      error[ii] = calc_error(bounds[ii], bounds[ii+1]);
    }
  }

  void set_verbose(bool v) { iverbose = v; }
  void set_partition_max_iterations(int max_it) {
    partition_max_iterations = max_it;
  }
  void set_line_search_max_iterations(int max_it) {
    line_search_max_iterations = max_it;
  }
  void set_partition_tolerance(ep_real pt) {
    partition_tolerance = pt;
  }
  void set_cubic_interpolation(ep_real ci) {
    cubic_interpolation = ci;
  }
  void set_resolution(ep_real res) {
    resolution = res;
  }

  void print_bounds(int ni, ep_real* bounds);
  void print_error(int ni, ep_real* error);


private:
  ep_real next_bound_above(ep_real target_error,
    ep_real bound1, ep_real boundn, 
    ep_real bound2_test, ep_real* error_test = 0);

  ep_real next_bound_below(ep_real target_error,
    ep_real bound0, ep_real bound2,
    ep_real bound1_test, ep_real* error_test = 0);

  EpStatus line_search(int ni, ep_real* bounds, ep_real* newbounds,
		       ep_real* error);

  ep_real next_bound_error_tolerance = 0.05;
  ep_real partition_tolerance = 0.05;
  ep_real resolution = 0.0; // Minimum meaningful distance between bounds
  int next_bound_max_iterations = 20;
  int partition_max_iterations = 20;
  int line_search_max_iterations = 10;

  bool iverbose = false;
  bool dump_error_iteration = false;
  bool cubic_interpolation = true;

};

#endif
