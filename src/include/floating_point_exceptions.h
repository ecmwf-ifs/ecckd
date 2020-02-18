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
