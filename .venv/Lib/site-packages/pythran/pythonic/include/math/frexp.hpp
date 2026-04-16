#ifndef PYTHONIC_INCLUDE_MATH_FREXP_HPP
#define PYTHONIC_INCLUDE_MATH_FREXP_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  std::tuple<double, long> frexp(double x);
  DEFINE_FUNCTOR(pythonic::math, frexp);
} // namespace math
PYTHONIC_NS_END

#endif
