#ifndef PYTHONIC_INCLUDE_NUMPY_FULL_HPP
#define PYTHONIC_INCLUDE_NUMPY_FULL_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class pS, class F, class dtype>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> full(pS const &shape, F fill_value,
                                                                 dtype d);

  template <class F, class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>> full(long size, F fill_value, dtype d);

  template <long N, class F, class dtype>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  full(std::integral_constant<long, N>, F fill_value, dtype d);

  template <class pS, class F>
  types::ndarray<F, sutils::shape_t<pS>> full(pS const &shape, F fill_value,
                                              types::none_type _ = {});

  template <class F>
  types::ndarray<F, types::pshape<long>> full(long size, F fill_value, types::none_type _ = {});

  template <long N, class F>
  types::ndarray<F, types::pshape<std::integral_constant<long, N>>>
  full(std::integral_constant<long, N>, F fill_value, types::none_type _ = {});

  DEFINE_FUNCTOR(pythonic::numpy, full);
} // namespace numpy
PYTHONIC_NS_END

#endif
