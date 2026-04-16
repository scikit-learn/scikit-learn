#ifndef PYTHONIC_NUMPY_ROLLAXIS_HPP
#define PYTHONIC_NUMPY_ROLLAXIS_HPP

#include "pythonic/include/numpy/rollaxis.hpp"

#include "pythonic/numpy/copy.hpp"
#include "pythonic/numpy/transpose.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
  rollaxis(types::ndarray<T, pS> const &a, long axis, long start)
  {
    long constexpr N = std::tuple_size<pS>::value;
    long t[N];
    if (start >= axis) {
      for (long i = 0; i < axis; ++i)
        t[i] = i;
      for (long i = axis + 1; i < start; ++i)
        t[i - 1] = i;
      t[start - 1] = axis;
      for (long i = start; i < N; ++i)
        t[i] = i;
    } else {
      for (long i = 0; i < start; ++i)
        t[i] = i;
      t[start] = axis;
      for (long i = start + 1; i <= axis; ++i)
        t[i] = i - 1;
      for (long i = axis + 1, n = N; i < n; ++i)
        t[i] = i;
    }
    return _transposer(a, t);
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(rollaxis);
} // namespace numpy
PYTHONIC_NS_END

#endif
