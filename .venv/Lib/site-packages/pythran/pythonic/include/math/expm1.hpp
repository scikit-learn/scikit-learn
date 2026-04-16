#ifndef PYTHONIC_INCLUDE_MATH_EXPM1_HPP
#define PYTHONIC_INCLUDE_MATH_EXPM1_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(expm1, std::expm1);
}
PYTHONIC_NS_END

#endif
