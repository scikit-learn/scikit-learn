#ifndef PYTHONIC_INCLUDE_NUMPY_EMPTY_HPP
#define PYTHONIC_INCLUDE_NUMPY_EMPTY_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  typename dtype::type empty(types::pshape<> const &shape, dtype d = dtype());

  template <class pS, class dtype = functor::float64>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> empty(pS const &shape,
                                                                  dtype d = dtype());

  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>> empty(long size, dtype d = dtype());

  template <long N, class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  empty(std::integral_constant<long, N>, dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, empty);
} // namespace numpy
PYTHONIC_NS_END

#endif
