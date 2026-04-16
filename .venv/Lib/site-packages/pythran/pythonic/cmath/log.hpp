#ifndef PYTHONIC_CMATH_LOG_HPP
#define PYTHONIC_CMATH_LOG_HPP

#include "pythonic/include/cmath/log.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  using std::log;
  double log(double x, double base)
  {
    return log(x) / log(base);
  }
} // namespace cmath
PYTHONIC_NS_END

#endif
