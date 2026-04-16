#ifndef PYTHONIC_INCLUDE_CMATH_E_HPP
#define PYTHONIC_INCLUDE_CMATH_E_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  double constexpr e = std::exp(1);
}
PYTHONIC_NS_END

#endif
