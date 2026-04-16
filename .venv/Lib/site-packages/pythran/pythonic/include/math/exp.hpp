#ifndef PYTHONIC_INCLUDE_MATH_EXP_HPP
#define PYTHONIC_INCLUDE_MATH_EXP_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(exp, std::exp);
}
PYTHONIC_NS_END

#endif
