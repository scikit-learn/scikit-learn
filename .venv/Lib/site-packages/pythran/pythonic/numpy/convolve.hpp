#ifndef PYTHONIC_NUMPY_CONVOLVE_HPP
#define PYTHONIC_NUMPY_CONVOLVE_HPP

#include "pythonic/include/numpy/convolve.hpp"
#include "pythonic/numpy/conjugate.hpp"
#include "pythonic/numpy/correlate.hpp"
#include "pythonic/numpy/flip.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class A, class B, typename U>
  types::ndarray<typename A::dtype, types::pshape<long>> convolve(A const &inA, B const &inB,
                                                                  U type)
  {
    auto inB_flipped = functor::flip{}(inB, 0);
    auto inB_flip_conj = functor::conjugate{}(inB_flipped);
    return functor::correlate{}(inA, inB_flip_conj, type);
  }

  template <class A, class B>
  types::ndarray<typename A::dtype, types::pshape<long>> convolve(A const &inA, B const &inB)
  {
    auto inB_flipped = functor::flip{}(inB, 0);
    auto inB_flip_conj = functor::conjugate{}(inB_flipped);
    return functor::correlate{}(inA, inB_flip_conj, "full");
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(convolve)
} // namespace numpy
PYTHONIC_NS_END

#endif
