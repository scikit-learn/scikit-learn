#ifndef PYTHONIC_NUMPY_ONES_HPP
#define PYTHONIC_NUMPY_ONES_HPP

#include "pythonic/include/numpy/ones.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype>
  typename dtype::type ones(std::tuple<> const &shape, dtype d)
  {
    return static_cast<typename dtype::type>(1);
  }

  template <class pS, class dtype>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> ones(pS const &shape, dtype d)
  {
    return {(sutils::shape_t<pS>)shape, typename dtype::type(1)};
  }

  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>> ones(long size, dtype d)
  {
    return ones(types::pshape<long>(size), d);
  }

  template <long N, class dtype>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  ones(std::integral_constant<long, N>, dtype d)
  {
    return ones(types::pshape<std::integral_constant<long, N>>({}), d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
