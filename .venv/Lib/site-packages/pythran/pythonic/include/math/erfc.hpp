#ifndef PYTHONIC_INCLUDE_MATH_ERFC_HPP
#define PYTHONIC_INCLUDE_MATH_ERFC_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(erfc, std::erfc);
}
PYTHONIC_NS_END

#endif
