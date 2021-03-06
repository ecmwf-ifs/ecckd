// constants.h - Define some physical constants
//
// Copyright (C) 2019- ECMWF.
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <adept_arrays.h>

static const std::string K_NAME            = "molar_absorption_coeff";
static const adept::Real ACCEL_GRAVITY     = 9.80665; // m s-2
static const adept::Real SPECIFIC_HEAT_AIR = 1004.0;  // J kg-1 K-1
static const adept::Real LW_DIFFUSIVITY    = 1.66;
static const adept::Real MOLAR_MASS_DRY_AIR= 28.970; // g mol-1
static const adept::Real AVOGADRO_CONSTANT = 6.02214076e23; // mol-1
#endif
