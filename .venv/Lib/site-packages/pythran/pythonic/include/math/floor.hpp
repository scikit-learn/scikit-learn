#ifndef PYTHONIC_INCLUDE_MATH_FLOOR_HPP
#define PYTHONIC_INCLUDE_MATH_FLOOR_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  long floor(T x);

  DEFINE_FUNCTOR(pythonic::math, floor);
} // namespace math
PYTHONIC_NS_END

#endif
