#ifndef PYTHONIC_MATH_FACTORIAL_HPP
#define PYTHONIC_MATH_FACTORIAL_HPP

#include "pythonic/include/math/factorial.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{

  template <class T>
  T factorial(T x)
  {
    long res = 1;
    for (long i = 2; i <= x; i++)
      res *= i;
    return res;
  }
} // namespace math
PYTHONIC_NS_END

#endif
