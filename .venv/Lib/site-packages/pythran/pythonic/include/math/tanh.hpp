#ifndef PYTHONIC_INCLUDE_MATH_TANH_HPP
#define PYTHONIC_INCLUDE_MATH_TANH_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(tanh, std::tanh);
}
PYTHONIC_NS_END

#endif
