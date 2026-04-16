#ifndef PYTHONIC_INCLUDE_MATH_ACOSH_HPP
#define PYTHONIC_INCLUDE_MATH_ACOSH_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(acosh, std::acosh);
}
PYTHONIC_NS_END

#endif
