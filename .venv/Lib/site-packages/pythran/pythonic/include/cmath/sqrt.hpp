#ifndef PYTHONIC_INCLUDE_CMATH_SQRT_HPP
#define PYTHONIC_INCLUDE_CMATH_SQRT_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(sqrt, std::sqrt);
}
PYTHONIC_NS_END

#endif
