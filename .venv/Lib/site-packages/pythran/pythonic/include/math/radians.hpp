#ifndef PYTHONIC_INCLUDE_MATH_RADIANS_HPP
#define PYTHONIC_INCLUDE_MATH_RADIANS_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/math/pi.hpp"

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  double radians(T x);
  DEFINE_FUNCTOR(pythonic::math, radians);
} // namespace math
PYTHONIC_NS_END

#endif
