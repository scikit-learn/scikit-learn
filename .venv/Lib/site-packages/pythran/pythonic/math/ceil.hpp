#ifndef PYTHONIC_MATH_CEIL_HPP
#define PYTHONIC_MATH_CEIL_HPP

#include "pythonic/include/math/ceil.hpp"

#include "pythonic/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  long ceil(T x)
  {
    return std::ceil(x);
  }
} // namespace math
PYTHONIC_NS_END

#endif
