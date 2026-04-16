#ifndef PYTHONIC_INCLUDE_CMATH_LOG10_HPP
#define PYTHONIC_INCLUDE_CMATH_LOG10_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(log10, std::log10);
}
PYTHONIC_NS_END

#endif
