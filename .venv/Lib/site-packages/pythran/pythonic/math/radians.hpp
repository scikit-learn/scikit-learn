#ifndef PYTHONIC_MATH_RADIANS_HPP
#define PYTHONIC_MATH_RADIANS_HPP

#include "pythonic/include/math/radians.hpp"

#include "pythonic/math/pi.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  double radians(T x)
  {
    return (x * 2. * pi) / 360.;
  }
} // namespace math
PYTHONIC_NS_END

#endif
