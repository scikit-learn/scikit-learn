#ifndef PYTHONIC_INCLUDE_NUMPY_CORRELATE_HPP
#define PYTHONIC_INCLUDE_NUMPY_CORRELATE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class A, class B>
  types::ndarray<typename A::dtype, types::pshape<long>>
  correlate(A const &inA, B const &inB, types::str const &renorm = types::str("valid"));

  NUMPY_EXPR_TO_NDARRAY0_DECL(correlate)
  DEFINE_FUNCTOR(pythonic::numpy, correlate)
} // namespace numpy
PYTHONIC_NS_END

#endif
