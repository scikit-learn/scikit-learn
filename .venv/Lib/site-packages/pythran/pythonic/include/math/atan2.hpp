#ifndef PYTHONIC_INCLUDE_MATH_ATAN2_HPP
#define PYTHONIC_INCLUDE_MATH_ATAN2_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(atan2, std::atan2);
}
PYTHONIC_NS_END

#endif
