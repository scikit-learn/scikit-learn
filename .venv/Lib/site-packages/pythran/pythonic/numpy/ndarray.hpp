#ifndef PYTHONIC_NUMPY_NDARRAY_HPP
#define PYTHONIC_NUMPY_NDARRAY_HPP

#include "pythonic/include/numpy/ndarray.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/nested_container.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace anonymous
  {

    template <class pS, class dtype>
    types::ndarray<typename dtype::type, sutils::shape_t<pS>> ndarray(pS const &shape, dtype)
    {
      return {(sutils::shape_t<pS>)shape, builtins::None};
    }

    template <class dtype>
    types::ndarray<typename dtype::type, types::pshape<long>> ndarray(long size, dtype d)
    {
      return ndarray(types::pshape<long>(size), d);
    }

    template <long N, class dtype>
    types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
    ndarray(std::integral_constant<long, N>, dtype d)
    {
      return ndarray(types::pshape<std::integral_constant<long, N>>({}), d);
    }

  } // namespace anonymous
} // namespace numpy
PYTHONIC_NS_END

#endif
