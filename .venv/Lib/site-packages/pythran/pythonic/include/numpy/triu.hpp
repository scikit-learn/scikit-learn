#ifndef PYTHONIC_INCLUDE_NUMPY_TRIU_HPP
#define PYTHONIC_INCLUDE_NUMPY_TRIU_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, pS> triu(types::ndarray<T, pS> const &expr, int k = 0);

  NUMPY_EXPR_TO_NDARRAY0_DECL(triu)
  DEFINE_FUNCTOR(pythonic::numpy, triu)
} // namespace numpy
PYTHONIC_NS_END

#endif
