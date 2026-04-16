#ifndef PYTHONIC_INCLUDE_CMATH_COS_HPP
#define PYTHONIC_INCLUDE_CMATH_COS_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  template <class T>
  std::complex<T> cos(std::complex<T> const &v);
  template <class T>
  std::complex<T> cos(T const &v);

  DEFINE_FUNCTOR(pythonic::cmath, cos);
} // namespace cmath
PYTHONIC_NS_END

#endif
