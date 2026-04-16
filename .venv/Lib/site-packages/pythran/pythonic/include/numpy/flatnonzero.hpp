#ifndef PYTHONIC_INCLUDE_NUMPY_FLATNONZERO_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLATNONZERO_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<long, types::pshape<long>> flatnonzero(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, flatnonzero);
} // namespace numpy
PYTHONIC_NS_END

#endif
