#ifndef PYTHONIC_INCLUDE_MATH_TAN_HPP
#define PYTHONIC_INCLUDE_MATH_TAN_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(tan, std::tan);
}
PYTHONIC_NS_END

#endif
