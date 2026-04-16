#ifndef PYTHONIC_INCLUDE_NUMPY_NDIM_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDIM_HPP

#include "pythonic/include/numpy/shape.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  long ndim(E const &e);

  DEFINE_FUNCTOR(pythonic::numpy, ndim)
} // namespace numpy
PYTHONIC_NS_END

#endif
