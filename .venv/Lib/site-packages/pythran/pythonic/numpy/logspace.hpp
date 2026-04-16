#ifndef PYTHONIC_NUMPY_LOGSPACE_HPP
#define PYTHONIC_NUMPY_LOGSPACE_HPP

#include "pythonic/include/numpy/logspace.hpp"

#include "pythonic/numpy/linspace.hpp"
#include "pythonic/numpy/power.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  auto logspace(double start, double stop, long num, bool endpoint, double base)
      -> decltype(functor::power()(base, functor::linspace()(start, stop, num, endpoint)))
  {
    return functor::power()(base, functor::linspace()(start, stop, num, endpoint));
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
