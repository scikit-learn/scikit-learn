#ifndef PYTHONIC_INCLUDE_MATH_LGAMMA_HPP
#define PYTHONIC_INCLUDE_MATH_LGAMMA_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(lgamma, std::lgamma);
}
PYTHONIC_NS_END

#endif
