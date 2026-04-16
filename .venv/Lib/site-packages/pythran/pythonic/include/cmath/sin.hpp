#ifndef PYTHONIC_INCLUDE_CMATH_SIN_HPP
#define PYTHONIC_INCLUDE_CMATH_SIN_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(sin, std::sin);
}
PYTHONIC_NS_END

#endif
