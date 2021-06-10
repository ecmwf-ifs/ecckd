// single_gas_data.h - Define SingleGasData used by find_g_points
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

// This file defines the SingleGasData type used by
// find_g_points_lw.cpp; unfortunately, ckd_model.h defines a
// different type of the same name.

#ifndef SINGLE_GAS_DATA_H
#define SINGLE_GAS_DATA_H 1

#include <vector>
#include <adept_arrays.h>

#include "Error.h"

using namespace adept;

// Structure for storing the properties of the g points for a single
// gas
struct SingleGasData {
  SingleGasData() {};
  SingleGasData(const std::string& molecule_,   // Name of molecule
		const intVector& n_g_points_,   // Number of g-points in each band
		const intVector& band_number_,  // Band to which each g-point belongs
		const intVector& rank1_,        // { Index of first and last sorted
		const intVector& rank2_,        // {   wavenumber in this g-point
		const Vector& error_,           // RMS heating-rate error for each g-point
		const Vector& sorting_variable_,// Median of variable used to sort in each g
		const intVector& rank_)         // Rank of each wavenumber
    // Use "eval" to ensure deep copy
    : molecule(molecule_), n_g_points(eval(n_g_points_)),
      band_number(eval(band_number_)), rank1(eval(rank1_)), rank2(eval(rank2_)),
      error(eval(error_)), sorting_variable(eval(sorting_variable_)) {
    Molecule = molecule;
    std::transform(Molecule.begin(), Molecule.end(), Molecule.begin(), ::toupper);
    store_g_points(rank_);
  }

  // Number of bands
  int nbands() const { return n_g_points.size(); }

  // Number of g-points
  int ng() const { return rank1.size(); }

  void store_g_points(const intVector& rank_) {
    g_point.resize(maxval(rank_)+1);
    g_point = -1;
    for (int ig = 0; ig < rank1.size(); ++ig) {
      g_point.where(rank_ >= rank1(ig) && rank_ <= rank2(ig)) = ig;
    }
  }

  void print() {
    LOG << "Single-gas data for " << Molecule << ":\n";
    LOG << "  number of g-points in each band     = " << n_g_points << "\n";
    LOG << "  band associated with each g-point   = " << band_number << "\n";
    LOG << "  wavenumber rank lower bound         = " << rank1 << "\n";
    LOG << "  wavenumber rank upper bound         = " << rank2 << "\n";
    LOG << "  heating-rate error for each g-point = " << error << "\n";
    LOG << "  sorting variable for each g-point   = " << sorting_variable << "\n";
  }

  void clear() {
    molecule.clear();
    Molecule.clear();
    n_g_points.clear();
    band_number.clear();
    rank1.clear();
    rank2.clear();
    g_min.clear();
    g_max.clear();
    error.clear();
    sorting_variable.clear();
    g_point.clear();
  }

  SingleGasData& operator=(const SingleGasData& in) {
    clear();
    molecule = in.molecule;
    Molecule = in.Molecule;
    n_g_points = in.n_g_points;
    band_number = in.band_number;
    rank1 = in.rank1;
    rank2 = in.rank2;
    g_min = in.g_min;
    g_max = in.g_max;
    error = in.error;
    sorting_variable = in.sorting_variable;
    g_point = in.g_point;
    return *this;
  }

  // Data
  std::string molecule, Molecule;
  intVector n_g_points, band_number;
  intVector rank1, rank2;
  intVector g_min, g_max;
  Vector error, sorting_variable;
  intVector g_point;
};

int overlap_g_points(std::vector<SingleGasData>& gas_data,
		     intVector& band_number);

// Repartition g-points so that they are more evenly spaced (often the
// most optically thick one covers too narrow a region of the spectrum
void repartition_g_points(const SingleGasData& src, // Initial g-point distribution
			  const Vector& weight,     // Weight of each wavenumber
			  const intVector& rank,    // Rank of each wavenumber
			  SingleGasData& dest,      // Destination g-point distribution
			  intVector n_g_points = intVector()); // Target number of g-points


#endif
