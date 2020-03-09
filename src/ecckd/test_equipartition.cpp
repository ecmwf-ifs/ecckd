#include <adept_arrays.h>
#include <cmath>
#include <iostream>

#include "equipartition.h"

using namespace adept;

class TestEquipartition : public Equipartition {
public:
  TestEquipartition(int n) : npoints(n) {
    values = exp(linspace(-2,10.0,npoints));
  }
  ep_real calc_error(ep_real bound1, ep_real bound2) {
    total_comp_cost += (bound2-bound1);
    int ibound1 = std::ceil(bound1*(npoints-1));
    int ibound2 = std::floor(bound2*(npoints-1));
    return abs(sum(values(range(ibound1,ibound2)))
	       - (ibound2-ibound1+1)*values((ibound1+ibound2)/2));
  }

  int npoints;
  Vector values;
  ep_real total_comp_cost;

};


int main()
{
  adept::set_array_print_style(PRINT_STYLE_PLAIN);
  ep_real e_save;

  bool iverbose = false;

  int npoints = 1000000;

  TestEquipartition te(npoints);
  te.set_verbose(iverbose);
  te.set_partition_max_iterations(200);
  te.set_line_search_max_iterations(15);
  te.set_partition_tolerance(0.001);
  te.set_resolution(1.0/npoints);
  te.total_comp_cost = 0.0;

  for (int ic = 0; ic < 2; ++ic)
  {

    te.set_cubic_interpolation(ic == 1);
    std::cout << "TESTING EQUIPARTITION SPECIFYING N\n";
    int ni = 16;
    Vector bounds(ni+1);
    Vector error(ni);

    bounds = linspace(0.0, 1.0, ni+1);
    EpStatus status = te.equipartition_n(ni, bounds.data(), error.data());
    std::cout << "***STATUS*** " << ep_status_string(status) << "\n";
    
    std::cout << "bounds = " << bounds << "\n";
    std::cout << "error  = " << error  << "\n";
    std::cout << "  computational cost = " << te.total_comp_cost << "\n";
    ep_real mean_error, cost_fn;
    ep_real frac_std, frac_range;
    ep_stats(ni, error.data(), mean_error, cost_fn,
	     &frac_std, &frac_range);
    std::cout << "  mean error = " << mean_error
	      << "\n  cost function = " << cost_fn
	      << "\n  frac std = " << frac_std
	      << "\n  frac range = " << frac_range << "\n";
    e_save = error[0];
  }

  std::cout << "\n";

  if (0)
  { 
    std::cout << "TESTING EQUIPARTITION SPECIFYING E\n";
    int ni_max = 100;
    std::vector<ep_real> bounds;
    std::vector<ep_real> error;
    int ni = ni_max;
    te.total_comp_cost = 0.0;
    EpStatus status = te.equipartition_e(e_save, 0.0, 1.0, ni,
					 bounds, error);
    std::cout << "***STATUS*** " << ep_status_string(status) << "\n";
    std::cout << "bounds = " << Vector(&bounds[0],dimensions(ni+1)) << "\n";
    std::cout << "error  = " << Vector(&error[0],dimensions(ni))  << "\n";
    std::cout << "  computational cost = " << te.total_comp_cost << "\n";
    ep_real mean_error, cost_fn;
    ep_real frac_std, frac_range;
    ep_stats(ni, &error[0], mean_error, cost_fn,
	     &frac_std, &frac_range);
    std::cout << "  mean error = " << mean_error
	      << "\n  cost function = " << cost_fn
	      << "\n  frac std = " << frac_std
	      << "\n  frac range = " << frac_range << "\n";
  }

}
