#ifndef PYTHONIC_INCLUDE_NUMPY_AVERAGE_HPP
#define PYTHONIC_INCLUDE_NUMPY_AVERAGE_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/numpy/sum.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto average(E const &expr, types::none_type const &axis = builtins::None)
      -> decltype(sum(expr, axis) / 1.);

  template <class E>
  auto average(E const &expr, long axis) -> decltype(sum(expr, axis) / 1.);

  template <class E, class W>
  auto average(E const &expr, types::none_type const &axis, W const &weights)
      -> decltype(average(expr * asarray(weights) / average(asarray(weights))));

  DEFINE_FUNCTOR(pythonic::numpy, average);
} // namespace numpy
PYTHONIC_NS_END

#endif
