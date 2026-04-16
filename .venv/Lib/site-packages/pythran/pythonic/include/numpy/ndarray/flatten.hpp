#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_FLATTEN_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_FLATTEN_HPP

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, class pS>
    types::ndarray<T, types::pshape<long>> flatten(types::ndarray<T, pS> const &a);

    NUMPY_EXPR_TO_NDARRAY0_DECL(flatten);
    DEFINE_FUNCTOR(pythonic::numpy::ndarray, flatten);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
