#ifndef PYTHONIC_NUMPY_ZEROS_HPP
#define PYTHONIC_NUMPY_ZEROS_HPP

#include "pythonic/include/numpy/zeros.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype>
  typename dtype::type zeros(std::tuple<> const &shape, dtype d)
  {
    return static_cast<typename dtype::type>(0);
  }

  template <class pS, class dtype>
  types::ndarray<typename dtype::type, sutils::shape_t<pS>> zeros(pS const &shape, dtype d)
  {
    using T = typename dtype::type;
    // use calloc even if we have a non integer type. This looks ok on modern
    // architecture, although not really standard
    auto *buffer = utils::callocate<T>(sutils::sprod(shape));
    return {buffer, (sutils::shape_t<pS>)shape, types::ownership::owned};
  }

  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>> zeros(long size, dtype d)
  {
    return zeros(types::pshape<long>(size), d);
  }

  template <long N, class dtype>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
  zeros(std::integral_constant<long, N>, dtype d)
  {
    return zeros(types::pshape<std::integral_constant<long, N>>({}), d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
