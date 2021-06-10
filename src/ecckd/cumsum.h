// cumsum.h - Matlab-like cumulative summation for Adept
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

#ifndef Cumsum_H
#define Cumsum_H 1

#include <adept_arrays.h>

// Cumulative summation for 1D arrays
template <typename T, bool IsActive>
inline
adept::Array<1,T,IsActive> cumsum(const adept::Array<1,T,IsActive> rhs) {
  adept::Array<1,T,IsActive> ans(rhs.size());
  if (!rhs.empty()) {
    ans(0) = rhs(0);
    for (adept::Index ii = 1; ii < rhs.size(); ++ii) {
      ans(ii) = ans(ii-1) + rhs(ii);
    }
  }
  return ans;
}


#endif
