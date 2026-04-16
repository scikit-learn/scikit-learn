#ifndef PYTHONIC_NUMPY_IDENTITY_HPP
#define PYTHONIC_NUMPY_IDENTITY_HPP

#include "pythonic/include/numpy/identity.hpp"

#include "pythonic/numpy/eye.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype>
  auto identity(long n, dtype d) -> decltype(eye(n, n, 0, d))
  {
    return eye(n, n, 0, d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
