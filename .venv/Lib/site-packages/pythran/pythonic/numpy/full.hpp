#ifndef PYTHONIC_NUMPY_FULL_HPP
#define PYTHONIC_NUMPY_FULL_HPP

#include "pythonic/include/numpy/full.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class pS, class F, class dtype>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> full(pS const &shape, F fill_value,
                                                                 dtype d)
  {
    return {(sutils::shape_t<pS>)shape, typename dtype::type(fill_value)};
  }

  template <class F, class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>> full(long size, F fill_value, dtype d)
  {
    return full(types::pshape<long>(size), fill_value, d);
  }

  template <long N, class F, class dtype>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  full(std::integral_constant<long, N>, F fill_value, dtype d)
  {
    return full(types::pshape<std::integral_constant<long, N>>({}), fill_value, d);
  }

  template <class pS, class F>
  types::ndarray<F, sutils::shape_t<pS>> full(pS const &shape, F fill_value, types::none_type)
  {
    return {(sutils::shape_t<pS>)shape, fill_value};
  }

  template <class F>
  types::ndarray<F, types::pshape<long>> full(long size, F fill_value, types::none_type nt)
  {
    return full(types::pshape<long>(size), fill_value, nt);
  }

  template <long N, class F>
  types::ndarray<F, types::pshape<std::integral_constant<long, N>>>
  full(std::integral_constant<long, N>, F fill_value, types::none_type nt)
  {
    return full(types::pshape<std::integral_constant<long, N>>({}), fill_value, nt);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
