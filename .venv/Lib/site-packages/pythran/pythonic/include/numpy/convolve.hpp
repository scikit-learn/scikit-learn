#ifndef PYTHONIC_INCLUDE_NUMPY_CONVOLVE_HPP
#define PYTHONIC_INCLUDE_NUMPY_CONVOLVE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class A, class B, typename U>
  types::ndarray<typename A::dtype, types::pshape<long>> convolve(A const &inA, B const &inB,
                                                                  U renorm = types::str("full"));

  template <class A, class B>
  types::ndarray<typename A::dtype, types::pshape<long>> convolve(A const &inA, B const &inB);

  NUMPY_EXPR_TO_NDARRAY0_DECL(convolve)
  DEFINE_FUNCTOR(pythonic::numpy, convolve)
} // namespace numpy
PYTHONIC_NS_END

#endif
