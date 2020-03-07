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
  EP_FAILURE,
  EP_INPUT_ERROR // Error in initial bounds (e.g. not monotonic)
} EpStatus;

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


private:
  ep_real next_bound_above(ep_real target_error,
    ep_real bound1, ep_real boundn, 
    ep_real bound2_test, ep_real* error_test = 0);
  ep_real next_bound_below(ep_real target_error,
    ep_real bound0, ep_real bound2,
    ep_real bound1_test, ep_real* error_test = 0);

  ep_real next_bound_error_tolerance = 0.05;
  ep_real partition_tolerance = 0.05;
  ep_real next_bound_max_iterations = 20;
  ep_real partition_max_iterations = 20;

};

#endif
