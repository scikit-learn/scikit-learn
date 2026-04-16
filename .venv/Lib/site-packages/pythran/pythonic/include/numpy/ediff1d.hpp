#ifndef PYTHONIC_INCLUDE_NUMPY_EDIFF1D_HPP
#define PYTHONIC_INCLUDE_NUMPY_EDIFF1D_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<typename E::dtype, types::pshape<long>> ediff1d(E const &expr);

  template <class E>
  auto ediff1d(types::list<E> const &expr) -> decltype(ediff1d(asarray(expr)));

  DEFINE_FUNCTOR(pythonic::numpy, ediff1d);
} // namespace numpy
PYTHONIC_NS_END

#endif
