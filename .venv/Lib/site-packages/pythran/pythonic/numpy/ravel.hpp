#ifndef PYTHONIC_NUMPY_RAVEL_HPP
#define PYTHONIC_NUMPY_RAVEL_HPP

#include "pythonic/include/numpy/ravel.hpp"

#include "pythonic/numpy/ndarray/reshape.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::pshape<long>> ravel(types::ndarray<T, pS> const &expr)
  {
    return expr.reshape(types::pshape<long>{expr.flat_size()});
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(ravel);
} // namespace numpy
PYTHONIC_NS_END

#endif
