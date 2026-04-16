#ifndef PYTHONIC_INCLUDE_NUMPY_RAVEL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RAVEL_HPP

#include "pythonic/include/numpy/ndarray/reshape.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::pshape<long>> ravel(types::ndarray<T, pS> const &expr);

  NUMPY_EXPR_TO_NDARRAY0_DECL(ravel);
  DEFINE_FUNCTOR(pythonic::numpy, ravel);
} // namespace numpy
PYTHONIC_NS_END

#endif
