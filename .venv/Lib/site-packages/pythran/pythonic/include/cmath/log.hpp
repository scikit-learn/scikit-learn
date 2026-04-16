#ifndef PYTHONIC_INCLUDE_CMATH_LOG_HPP
#define PYTHONIC_INCLUDE_CMATH_LOG_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  using std::log;
  double log(double x, double base);
  DEFINE_FUNCTOR(pythonic::cmath, log);
} // namespace cmath
PYTHONIC_NS_END

#endif
