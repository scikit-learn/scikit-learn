#ifndef PYTHONIC_INCLUDE_NUMPY_IDENTITY_HPP
#define PYTHONIC_INCLUDE_NUMPY_IDENTITY_HPP

#include "pythonic/include/numpy/eye.hpp"
#include "pythonic/include/numpy/float64.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype = functor::float64>
  auto identity(long n, dtype d = dtype()) -> decltype(eye(n, n, 0, d));

  DEFINE_FUNCTOR(pythonic::numpy, identity);
} // namespace numpy
PYTHONIC_NS_END

#endif
