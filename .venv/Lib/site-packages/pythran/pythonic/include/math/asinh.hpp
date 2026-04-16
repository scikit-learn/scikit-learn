#ifndef PYTHONIC_INCLUDE_MATH_ASINH_HPP
#define PYTHONIC_INCLUDE_MATH_ASINH_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(asinh, std::asinh);
}
PYTHONIC_NS_END

#endif
