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
