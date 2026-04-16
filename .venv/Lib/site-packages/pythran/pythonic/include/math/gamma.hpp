#ifndef PYTHONIC_INCLUDE_MATH_GAMMA_HPP
#define PYTHONIC_INCLUDE_MATH_GAMMA_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  double gamma(double x);
  DEFINE_FUNCTOR(pythonic::math, gamma);
} // namespace math
PYTHONIC_NS_END

#endif
