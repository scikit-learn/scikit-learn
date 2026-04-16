#ifndef PYTHONIC_INCLUDE_MATH_LOG1P_HPP
#define PYTHONIC_INCLUDE_MATH_LOG1P_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(log1p, std::log1p);
}
PYTHONIC_NS_END

#endif
