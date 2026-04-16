#ifndef PYTHONIC_INCLUDE_MATH_FACTORIAL_HPP
#define PYTHONIC_INCLUDE_MATH_FACTORIAL_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{

  template <class T>
  T factorial(T x);

  DEFINE_FUNCTOR(pythonic::math, factorial);
} // namespace math
PYTHONIC_NS_END

#endif
