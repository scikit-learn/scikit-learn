#ifndef PYTHONIC_INCLUDE_CMATH_PI_HPP
#define PYTHONIC_INCLUDE_CMATH_PI_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  double constexpr pi = std::atan(1) * 4;
}
PYTHONIC_NS_END

#endif
