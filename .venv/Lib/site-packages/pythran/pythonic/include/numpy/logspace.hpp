#ifndef PYTHONIC_INCLUDE_NUMPY_LOGSPACE_HPP
#define PYTHONIC_INCLUDE_NUMPY_LOGSPACE_HPP

#include "pythonic/include/numpy/linspace.hpp"
#include "pythonic/include/numpy/power.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  auto logspace(double start, double stop, long num = 50, bool endpoint = true, double base = 10.0)
      -> decltype(functor::power()(base, functor::linspace()(start, stop, num, endpoint)));

  DEFINE_FUNCTOR(pythonic::numpy, logspace);
} // namespace numpy
PYTHONIC_NS_END

#endif
