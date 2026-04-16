#ifndef PYTHONIC_INCLUDE_MATH_SIN_HPP
#define PYTHONIC_INCLUDE_MATH_SIN_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(sin, std::sin);
}
PYTHONIC_NS_END

#endif
