#ifndef PYTHONIC_INCLUDE_MATH_POW_HPP
#define PYTHONIC_INCLUDE_MATH_POW_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(pow, std::pow);
}
PYTHONIC_NS_END

#endif
