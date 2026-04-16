#ifndef PYTHONIC_INCLUDE_MATH_ACOS_HPP
#define PYTHONIC_INCLUDE_MATH_ACOS_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(acos, std::acos);
}
PYTHONIC_NS_END

#endif
