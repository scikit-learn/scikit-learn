#ifndef PYTHONIC_NUMPY_EYE_HPP
#define PYTHONIC_NUMPY_EYE_HPP

#include "pythonic/include/numpy/eye.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/zeros.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype>
  types::ndarray<typename dtype::type, types::array_tuple<long, 2>> eye(long N, long M, long k,
                                                                        dtype d)
  {
    types::ndarray<typename dtype::type, types::array_tuple<long, 2>> out =
        zeros(types::make_tuple(N, M), d);
    if (k >= 0)
      for (int i = 0, j = k; i < N && j < M; ++i, ++j)
        out[i][j] = typename dtype::type(1);
    else
      for (int i = -k, j = 0; i < N && j < M; ++i, ++j)
        out[i][j] = typename dtype::type(1);
    return out;
  }

  template <class dtype>
  types::ndarray<typename dtype::type, types::array_tuple<long, 2>> eye(long N, types::none_type M,
                                                                        long k, dtype d)
  {
    return eye(N, N, k, d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
