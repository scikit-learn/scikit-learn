#ifndef PYTHONIC_NUMPY_TAKE_HPP
#define PYTHONIC_NUMPY_TAKE_HPP

#include "pythonic/include/numpy/take.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class F, class T>
  auto take(T &&expr, F &&indices) -> decltype(std::forward<T>(expr)[std::forward<T>(indices)])
  {
    return expr[indices];
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
