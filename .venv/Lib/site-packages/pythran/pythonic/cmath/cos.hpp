#ifndef PYTHONIC_CMATH_COS_HPP
#define PYTHONIC_CMATH_COS_HPP

#include "pythonic/include/cmath/cos.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  template <class T>
  std::complex<T> cos(std::complex<T> const &v)
  {
    return std::cos(v);
  }

  template <class T>
  std::complex<T> cos(T const &v)
  {
    return std::cos(v);
  }
} // namespace cmath
PYTHONIC_NS_END

#endif
