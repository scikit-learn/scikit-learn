#ifndef PYTHONIC_INCLUDE_MATH_ASIN_HPP
#define PYTHONIC_INCLUDE_MATH_ASIN_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(asin, std::asin);
}
PYTHONIC_NS_END

#endif
