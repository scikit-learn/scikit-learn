#ifndef PYTHONIC_INCLUDE_MATH_LOG_HPP
#define PYTHONIC_INCLUDE_MATH_LOG_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  using std::log;
  double log(double x, double base);
  DEFINE_FUNCTOR(pythonic::math, log);
} // namespace math
PYTHONIC_NS_END

#endif
