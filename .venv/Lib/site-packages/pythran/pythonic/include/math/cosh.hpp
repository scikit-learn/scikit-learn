#ifndef PYTHONIC_INCLUDE_MATH_COSH_HPP
#define PYTHONIC_INCLUDE_MATH_COSH_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(cosh, std::cosh);
}
PYTHONIC_NS_END

#endif
