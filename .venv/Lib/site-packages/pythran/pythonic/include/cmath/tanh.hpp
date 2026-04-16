#ifndef PYTHONIC_INCLUDE_CMATH_TANH_HPP
#define PYTHONIC_INCLUDE_CMATH_TANH_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(tanh, std::tanh);
}
PYTHONIC_NS_END

#endif
