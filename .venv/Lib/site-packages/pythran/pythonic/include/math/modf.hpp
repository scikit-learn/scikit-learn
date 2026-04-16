#ifndef PYTHONIC_INCLUDE_MATH_MODF_HPP
#define PYTHONIC_INCLUDE_MATH_MODF_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  std::tuple<double, double> modf(double x);
  DEFINE_FUNCTOR(pythonic::math, modf);
} // namespace math
PYTHONIC_NS_END

#endif
