/// @file      arrays.h
/// @brief     Includes the Adept header and implements some useful vector/matrix functions
/// @author    Robin J. Hogan
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef Arrays_H
#define Arrays_H 1

#include <adept_arrays.h>
#include "Error.h"

using namespace adept;

/// Return a vector containing the common elements, assuming both v1
/// and v2 are ordered
inline
adept::IntVector
common_elements(const adept::IntVector& v1, const adept::IntVector& v2)
{
  using namespace adept;
  IntVector v(v1.size());
  int i = 0, i2 = 0;
  for (int i1 = 0; i1 < v1.size() && i2 < v2.size(); ++i1) {
    while (i2 < v2.size()-1 && v2[i2] < v1[i1]) {
      ++i2;
    }
    if (v1[i1] == v2[i2]) {
      v[i++] = v2[i2++];
    }
  }
  if (i == 0) {
    return IntVector();
  }
  else {
    return v(range(0, i-1));
  }
}

/*
inline
adept::Vector
linspace(adept::Real start, adept::Real end, int n)
{
  using namespace adept;
  Vector out(n);
  Real d = (end-start)/(n-1);
  for (int j = 0; j < n; j++) {
    out(j) = start + d*j;
  }
  return out;
}
*/

inline
adept::Vector
logspace(adept::Real start, adept::Real end, int n)
{
  using namespace adept;
  Vector out(n);
  Real log_start = log(start);
  Real d = (log(end)-log_start)/(n-1);
  for (int j = 0; j < n; j++) {
    out(j) = exp(log_start + d*j);
  }
  return out;
}

/// Return an index to the first "true" element of boolean expression
/// x, or -1 if there are none
inline
int
find_first(const boolVector& x)
{
  for (int j = 0; j < x.size(); j++) {
    if (x(j)) {
      return j;
    }
  }
  return -1;
}

/// Return an index to the last "true" element of boolean expression
/// x, or -1 if there are none
inline
int
find_last(const boolVector& x)
{
  for (int j = x.size()-1; j >= 1; j--) {
    if (x(j)) {
      return j;
    }
  }
  return -1;
}

/*
inline
Vector
solve_least_squares(const Matrix& A, const Vector& b)
{
  if (A.dimension(0) != b.size()) {
    ERROR << "Attempt to perform least squares when number of rows in matrix A["
      << A.dimension(0) << "," << A.dimension(1) << "] and vector b[" << b.size() << "] do not match\n";
    THROW(UNEXPECTED_EXCEPTION);
  }
  int nx = A.dimension(1);
  Vector rhs(nx);
  SymmetricMatrix X(nx,nx);

  rhs = A.T() ** b;
  X = A.T() ** A;
  return solve(X, rhs);
}
*/

#endif
