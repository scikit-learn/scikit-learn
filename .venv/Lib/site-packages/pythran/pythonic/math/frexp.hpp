#ifndef PYTHONIC_MATH_FREXP_HPP
#define PYTHONIC_MATH_FREXP_HPP

#include "pythonic/include/math/frexp.hpp"

#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  std::tuple<double, long> frexp(double x)
  {
    int exp;
    double sig = std::frexp(x, &exp);
    return std::tuple<double, long>(sig, exp);
  }
} // namespace math
PYTHONIC_NS_END

#endif
