#ifndef PYTHONIC_INCLUDE_NUMPY_TRI_HPP
#define PYTHONIC_INCLUDE_NUMPY_TRI_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long, long>>
  tri(long N, long M = -1, long k = 0, dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, tri)
} // namespace numpy
PYTHONIC_NS_END

#endif
