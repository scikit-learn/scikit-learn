#ifndef PYTHONIC_INCLUDE_MATH_SINH_HPP
#define PYTHONIC_INCLUDE_MATH_SINH_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(sinh, std::sinh);
}
PYTHONIC_NS_END

#endif
