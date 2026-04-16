#ifndef PYTHONIC_NUMPY_TRI_HPP
#define PYTHONIC_NUMPY_TRI_HPP

#include "pythonic/include/numpy/tri.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long, long>> tri(long N, long M, long k,
                                                                      dtype d)
  {
    if (M == -1)
      M = N;
    types::ndarray<typename dtype::type, types::pshape<long, long>> out(
        types::pshape<long, long>{N, M}, 0);
    for (int i = 0; i < N; ++i)
      for (long j = 0; j < M; ++j)
        if (j - i <= k)
          out[i][j] = 1;
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
