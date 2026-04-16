#ifndef PYTHONIC_INCLUDE_MATH_ISNAN_HPP
#define PYTHONIC_INCLUDE_MATH_ISNAN_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(isnan, std::isnan);
}
PYTHONIC_NS_END

#endif
