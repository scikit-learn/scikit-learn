#ifndef PYTHONIC_NUMPY_TRACE_HPP
#define PYTHONIC_NUMPY_TRACE_HPP

#include "pythonic/include/numpy/trace.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T>
  typename T::dtype trace(T const &expr, int offset)
  {
    static_assert(T::value == 2, "Not Implemented : Trace for dimension != 2");

    typename T::dtype res = 0;
    long y_offset = std::max(-offset, 0);
    long x_offset = std::max(0, offset);
    long size = std::min(expr.flat_size() - y_offset, expr.fast(0).flat_size() - x_offset);
    if (offset < 0)
      for (long i = 0; i < size; i++)
        res += expr.fast(i + offset).fast(i);
    else if (offset > 0)
      for (long i = 0; i < size; i++)
        res += expr.fast(i).fast(i + offset);
    else
      for (long i = 0; i < size; i++)
        res += expr.fast(i).fast(i);
    return res;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
