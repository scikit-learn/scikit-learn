#ifndef PYTHONIC_MATH_TRUNC_HPP
#define PYTHONIC_MATH_TRUNC_HPP

#include "pythonic/include/math/trunc.hpp"

#include "pythonic/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  long trunc(T x)
  {
    return x;
  }
} // namespace math
PYTHONIC_NS_END

#endif
