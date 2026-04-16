#ifndef PYTHONIC_INCLUDE_MATH_CEIL_HPP
#define PYTHONIC_INCLUDE_MATH_CEIL_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  long ceil(T x);

  DEFINE_FUNCTOR(pythonic::math, ceil);
} // namespace math
PYTHONIC_NS_END

#endif
