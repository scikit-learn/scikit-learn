#ifndef PYTHONIC_MATH_MODF_HPP
#define PYTHONIC_MATH_MODF_HPP

#include "pythonic/include/math/modf.hpp"

#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{

  std::tuple<double, double> modf(double x)
  {
    double i;
    double frac = std::modf(x, &i);
    return std::make_tuple(frac, i);
  }
} // namespace math
PYTHONIC_NS_END

#endif
