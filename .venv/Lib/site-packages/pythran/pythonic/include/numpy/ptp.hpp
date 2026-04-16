#ifndef PYTHONIC_INCLUDE_NUMPY_PTP_HPP
#define PYTHONIC_INCLUDE_NUMPY_PTP_HPP

#include "pythonic/include/numpy/max.hpp"
#include "pythonic/include/numpy/min.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto ptp(E const &expr, long axis) -> decltype(max(expr, axis) - min(expr, axis));

  template <class E>
  auto ptp(E const &expr) -> decltype(max(expr) - min(expr));

  DEFINE_FUNCTOR(pythonic::numpy, ptp);
} // namespace numpy
PYTHONIC_NS_END

#endif
