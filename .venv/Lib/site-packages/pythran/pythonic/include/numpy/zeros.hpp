#ifndef PYTHONIC_INCLUDE_NUMPY_ZEROS_HPP
#define PYTHONIC_INCLUDE_NUMPY_ZEROS_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  typename dtype::type zeros(std::tuple<> const &shape, dtype d = dtype());

  template <class pS, class dtype = functor::float64>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> zeros(pS const &shape,
                                                                  dtype d = dtype());

  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>> zeros(long size, dtype d = dtype());

  template <long N, class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  zeros(std::integral_constant<long, N>, dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, zeros);
} // namespace numpy
PYTHONIC_NS_END

#endif
