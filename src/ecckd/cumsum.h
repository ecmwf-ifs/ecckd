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
