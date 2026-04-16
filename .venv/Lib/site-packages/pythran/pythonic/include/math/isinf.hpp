#ifndef PYTHONIC_INCLUDE_MATH_ISINF_HPP
#define PYTHONIC_INCLUDE_MATH_ISINF_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  bool isinf(T const &v)
  {
    return std::isinf(v);
  }
  DEFINE_FUNCTOR(pythonic::math, isinf);
} // namespace math
PYTHONIC_NS_END

#endif
