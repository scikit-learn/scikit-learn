#ifndef PYTHONIC_INCLUDE_MATH_LDEXP_HPP
#define PYTHONIC_INCLUDE_MATH_LDEXP_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(ldexp, std::ldexp);
}
PYTHONIC_NS_END

#endif
