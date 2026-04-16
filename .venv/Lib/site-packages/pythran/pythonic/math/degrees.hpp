#ifndef PYTHONIC_MATH_DEGREES_HPP
#define PYTHONIC_MATH_DEGREES_HPP

#include "pythonic/include/math/degrees.hpp"

#include "pythonic/math/pi.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{

  template <class T>
  double degrees(T x)
  {
    return (x * 360.) / (2. * pi);
  }
} // namespace math
PYTHONIC_NS_END

#endif
