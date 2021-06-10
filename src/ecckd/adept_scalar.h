// adept_scalar.h - Select an active or passive Adept scalar
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

#ifndef ADEPT_SCALAR_H
#define ADEPT_SCALAR_H 1

#include <adept_arrays.h>

/// Select an active versus static scalar depending on bool template
/// argument
template <bool IsActive>
struct scalar {
  typedef adept::Real type;
};

template <>
struct scalar<true> {
  typedef adept::aReal type;
};

#endif
