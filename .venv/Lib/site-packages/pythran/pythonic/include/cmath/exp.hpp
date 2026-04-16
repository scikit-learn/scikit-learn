#ifndef PYTHONIC_INCLUDE_CMATH_EXP_HPP
#define PYTHONIC_INCLUDE_CMATH_EXP_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(exp, std::exp);
}
PYTHONIC_NS_END

#endif
