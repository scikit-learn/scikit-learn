#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/nested_container.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace anonymous
  {

    template <class pS, class dtype = functor::float64>
    types::ndarray<typename dtype::type, sutils::shape_t<pS>> ndarray(pS const &shape,
                                                                      dtype d = dtype());

    template <class dtype = functor::float64>
    types::ndarray<typename dtype::type, types::pshape<long>> ndarray(long size, dtype d = dtype());

    template <long N, class dtype = functor::float64>
    types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, N>>>
    ndarray(std::integral_constant<long, N>, dtype d = dtype());

  } // namespace anonymous

  DEFINE_FUNCTOR(pythonic::numpy::anonymous, ndarray);
} // namespace numpy
PYTHONIC_NS_END

#endif
