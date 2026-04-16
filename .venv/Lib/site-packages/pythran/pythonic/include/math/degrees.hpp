#ifndef PYTHONIC_INCLUDE_MATH_DEGREES_HPP
#define PYTHONIC_INCLUDE_MATH_DEGREES_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{

  template <class T>
  double degrees(T x);

  DEFINE_FUNCTOR(pythonic::math, degrees);
} // namespace math
PYTHONIC_NS_END

#endif
