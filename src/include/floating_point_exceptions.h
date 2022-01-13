// floating_point_exceptions.h - Enable/disable floating-point trapping
//
// Copyright (C) 2019- ECMWF
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

#ifndef FLOATING_POINT_EXCEPTIONS
#define FLOATING_POINT_EXCEPTIONS 1

#include <fenv.h>

inline void enable_floating_point_exceptions() {
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
inline void disable_floating_point_exceptions() {
  feclearexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

#endif
