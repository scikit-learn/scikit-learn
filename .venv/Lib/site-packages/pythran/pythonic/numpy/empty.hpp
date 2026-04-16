#ifndef PYTHONIC_NUMPY_EMPTY_HPP
#define PYTHONIC_NUMPY_EMPTY_HPP

#include "pythonic/include/numpy/empty.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype>
  typename dtype::type empty(types::pshape<> const &shape, dtype)
  {
    return {};
  }

  template <class pS, class dtype>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> empty(pS const &shape, dtype)
  {
    return {(sutils::shape_t<pS>)shape, builtins::None};
  }

  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>> empty(long size, dtype d)
  {
    return empty(types::pshape<long>(size), d);
  }

  template <long N, class dtype>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  empty(std::integral_constant<long, N>, dtype d)
  {
    return empty(types::pshape<std::integral_constant<long, N>>({}), d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
