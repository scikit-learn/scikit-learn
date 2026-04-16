#ifndef PYTHONIC_INCLUDE_MATH_COS_HPP
#define PYTHONIC_INCLUDE_MATH_COS_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(cos, std::cos);
}
PYTHONIC_NS_END

#endif
