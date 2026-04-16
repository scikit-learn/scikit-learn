#ifndef PYTHONIC_INCLUDE_CMATH_ISNAN_HPP
#define PYTHONIC_INCLUDE_CMATH_ISNAN_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(isnan, std::isnan);
}
PYTHONIC_NS_END

#endif
