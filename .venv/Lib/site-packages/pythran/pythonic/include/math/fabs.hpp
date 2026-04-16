#ifndef PYTHONIC_INCLUDE_MATH_FABS_HPP
#define PYTHONIC_INCLUDE_MATH_FABS_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(fabs, std::fabs);
}
PYTHONIC_NS_END

#endif
