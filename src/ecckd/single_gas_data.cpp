// single_gas_data.cpp - Functions for manipulating g points
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

#include <cmath>

#include "single_gas_data.h"
#include "Error.h"

// Overlap the g-points of the various gases using the hypercube
// partition method of Hogan (2010)
int
overlap_g_points(std::vector<SingleGasData>& gas_data,
		 intVector& band_number) {
  const int ngas  = gas_data.size();
  const int nband = gas_data[0].n_g_points.size();
  intVector ng_band(nband);

  for (int iband = 0; iband < nband; ++iband) {
    // Number of g points per band is given by Eq. 7 of Hogan (2010):
    // 1-ngas+\sum_i=1^ngas(ng_i)
    ng_band(iband) = 1-ngas;
    for (int igas = 0; igas < ngas; ++igas) {
      ng_band(iband) += gas_data[igas].n_g_points[iband];
    }
  }

  int ng = sum(ng_band); // Total number of g points

  band_number.resize(ng);
  {
    int ig = 0;
    for (int iband = 0; iband < nband; ++iband) {
      band_number(range(ig,ig+ng_band(iband)-1)) = iband;
      ig += ng_band(iband);
    }
  }

  for (int igas = 0; igas < ngas; ++igas) {
    gas_data[igas].g_min.resize(ng);
    gas_data[igas].g_max.resize(ng);
  }

  int ig = 0; // Index of current "total" g point
  intVector ig_gas(ngas); // Index of current g point for each
			  // individual gas
  ig_gas = 0.0;

  for (int iband = 0; iband < nband; ++iband) {
    LOG << "Band " << iband << "\n";
    intVector ig_gas_start;
    ig_gas_start = ig_gas;
    // First "total" g point in each band is the intersection of the
    // first g point for each of the individual gases
    LOG << "  G-point " << ig << ": intersection of weakest spectral interval of each gas\n";
    for (int igas = 0; igas < ngas; ++igas) {
      gas_data[igas].g_min(ig) = ig_gas_start(igas);
      gas_data[igas].g_max(ig) = ig_gas_start(igas);
    }

    // Loop over total g points for this band
    for (int ig_band = 1; ig_band < ng_band(iband); ++ig_band) {
      // Find minimum sorting variable amongst the gases
      Real min_sorting_var = 1.0e30;
      int i_found_gas = -1;
      for (int igas = 0; igas < ngas; ++igas) {
	Real my_sorting_var = 1.0e30;
	if (ig_gas(igas) < ig_gas_start(igas)+gas_data[igas].n_g_points[iband]-1) {
	  my_sorting_var = gas_data[igas].sorting_variable(ig_gas(igas)+1);
	}
	//	else {
	//	  LOG << "  Reached maximum for gas " << igas << "\n";
	//	}
	if (my_sorting_var < min_sorting_var) {
	  min_sorting_var = my_sorting_var;
	  i_found_gas = igas;
	}
      }
      if (i_found_gas < 0) {
	ERROR << "Could not locate next gas to advance";
	throw;
      }

      ++ig_gas(i_found_gas);

      ++ig;

      LOG << "  G-point " << ig << ": major gas " << gas_data[i_found_gas].Molecule
	  << " (" << ig_gas(i_found_gas) << ")\n";

      for (int igas = 0; igas < ngas; ++igas) {
	if (i_found_gas == igas) {
	  gas_data[igas].g_min(ig) = ig_gas(igas);
	  gas_data[igas].g_max(ig) = ig_gas(igas);
	}
	else {
	  gas_data[igas].g_min(ig) = ig_gas_start(igas);
	  gas_data[igas].g_max(ig) = ig_gas(igas);
  	}
      } 
    }

    // Need to make sure that the indices are incremented ready for
    // the first "intersection of weakest" merged-gpoint of the next
    // band
    ++ig;
    for (int igas = 0; igas < ngas; ++igas) {
      ++ig_gas(igas);
    }

  }
  return ng;
}


// Repartition g-points so that they are more evenly spaced (often the
// most optically thick one covers too narrow a region of the spectrum
void
repartition_g_points(const SingleGasData& src, // Initial g-point distribution
		     const Vector& weight,     // Weight of each wavenumber
		     const intVector& rank,    // Rank of each wavenumber
		     SingleGasData& dest,      // Destination g-point distribution
		     intVector n_g_points)     // Target number of g-points
{
  dest.clear();

  dest.molecule = src.molecule;
  dest.Molecule = src.Molecule;

  if (n_g_points.empty()) {
    dest.n_g_points = src.n_g_points;
  }
  else {
    dest.n_g_points = n_g_points;
  }
  int nband = src.nbands();
  int ng    = sum(dest.n_g_points);

  // Work out band number corresponding to each new g-point
  dest.band_number.resize(ng);
  dest.rank1.resize(ng);
  dest.rank2.resize(ng);

  int igstart = 0;
  for (int iband = 0; iband < nband; ++iband) {
    dest.band_number(range(igstart,igstart-1+dest.n_g_points(iband))) = iband;
    igstart += dest.n_g_points(iband);
  }

  Vector weight_sorted(rank.size());
  weight_sorted(rank) = weight;

  Vector cum_error_density(rank.size());
  cum_error_density = 0.0;

  // Loop through each band, reordering g-points in each one
  int ioldg = 0;
  int ig = 0;

  for (int iband = 0; iband < nband; ++iband) {
    // Work out error density in each existing g-point
    Vector error_density(src.n_g_points(iband)); // mean error density
    Vector sum_weight(src.n_g_points(iband));
    if (iband > 0) {
      ioldg = sum(src.n_g_points(range(0,iband-1)));
    }
    dest.rank1(ig) = src.rank1(ioldg);
    for (int ioldglocal = 0; ioldglocal < src.n_g_points(iband); ++ioldglocal) {
      intVector index = find(rank >= src.rank1(ioldg) && rank <= src.rank2(ioldg));
      sum_weight(ioldglocal) = sum(weight(index));
      error_density(ioldglocal) = src.error(ioldg) / sum_weight(ioldglocal);
      ++ioldg;
    }

    // Error density at each end of g-point range
    Vector error_density1(src.n_g_points(iband));
    Vector error_density2(src.n_g_points(iband));
    // Flat for g-points at each end
    //    error_density1(0)   = error_density(0);
    //    error_density2(0)   = error_density(0);
    error_density1(end) = error_density(end);
    error_density2(end) = error_density(end);

    for (int ioldglocal = 0; ioldglocal < src.n_g_points(iband)-1; ++ioldglocal) {
      Real ideal_error_density1 = 0.0;
      if (ioldglocal > 0) {
	ideal_error_density1 = 0.5*(error_density(ioldglocal)+error_density(ioldglocal-1));
      }
      Real ideal_error_density2 = 0.5*(error_density(ioldglocal)+error_density(ioldglocal+1));
      if (ioldglocal == src.n_g_points(iband)-1) {
	ideal_error_density2 = error_density(end);
      }
      if ((ideal_error_density1 < error_density(ioldglocal))
	  == (error_density(ioldglocal) < ideal_error_density2)) {
	// error_density lies between ideal values
	Real diff = std::copysign(std::fmin(std::fabs(error_density(ioldglocal)-ideal_error_density1),
					    std::fabs(ideal_error_density2-error_density(ioldglocal))),
				  error_density(ioldglocal)-ideal_error_density1);
	error_density1(ioldglocal) = error_density(ioldglocal)-diff;
	error_density2(ioldglocal) = error_density(ioldglocal)+diff;
      }
      else {
	// error density does not lie between ideal values
	error_density1(ioldglocal) = error_density(ioldglocal);
	error_density2(ioldglocal) = error_density(ioldglocal);
      }
    }
    /*
    error_density1 = error_density;
    error_density2 = error_density;
    */

    LOG << "error_density = " << error_density << "\n";
    LOG << "error_density1 = " << error_density1 << "\n";
    LOG << "error_density2 = " << error_density2 << "\n";
    LOG << "sum_weight = " << sum_weight << "\n";

    LOG << "src.rank1 = " << src.rank1 << "\n";
    LOG << "src.rank2 = " << src.rank2 << "\n";

    LOG << "ng = " << ng << "\n";

    // Construct cumulative distribution of error
    Real sum_error_density = 0.0;
    if (iband > 0) {
      ioldg = sum(src.n_g_points(range(0,iband-1)));
    }
    else {
      ioldg = 0;
    }
    for (int ioldglocal = 0; ioldglocal < src.n_g_points(iband); ++ioldglocal, ++ioldg) {
      Real x = 0.0; // Cumulative normalized weight so far along old g point
      for (int irank = src.rank1(ioldg); irank <= src.rank2(ioldg); ++irank) {
	x += weight_sorted(irank)/sum_weight(ioldglocal);
	Real local_error_density = (1.0-x) * error_density1(ioldglocal)
	                               +x  * error_density2(ioldglocal);
	sum_error_density += weight_sorted(irank) * local_error_density;
	cum_error_density(irank) = sum_error_density;
      }
    }

    LOG << "sum_error_density = " << sum_error_density << "\n";

    // Repartition
    int irank = dest.rank1(ig);
    int iglocal = 0;
    while (iglocal < dest.n_g_points(iband)-1) {
      Real target = static_cast<Real>(iglocal+1)*sum_error_density
	/ static_cast<Real>(dest.n_g_points(iband));
      if (src.n_g_points(iband) == dest.n_g_points(iband)) {
	Real damper = 0.8;
	LOG << "FACTORS = " << sum(src.error(range(0,iglocal))) << " " << static_cast<Real>(iglocal+1)*sum_error_density
	  / static_cast<Real>(dest.n_g_points(iband)) << "\n";

	target = damper * sum(src.error(range(0,iglocal)))
	  + (1.0-damper) * target;
      }
      while (cum_error_density(irank) < target) {
	++irank;
      }
      dest.rank2(ig) = irank-1;
      ++ig;
      ++iglocal;
      dest.rank1(ig) = irank;
    }
    dest.rank2(ig) = src.rank2(sum(src.n_g_points(range(0,iband)))-1);
  }

  dest.store_g_points(rank);
  dest.error.resize(ng); dest.error = -1.0;
  dest.sorting_variable.resize(ng); dest.sorting_variable = -1.0;

}
