#ifndef PYTHONIC_INCLUDE_CMATH_TAN_HPP
#define PYTHONIC_INCLUDE_CMATH_TAN_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(tan, std::tan);
}
PYTHONIC_NS_END

#endif
