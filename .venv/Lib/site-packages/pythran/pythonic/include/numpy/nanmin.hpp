#ifndef PYTHONIC_INCLUDE_NUMPY_NANMIN_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANMIN_HPP

#include "pythonic/include/numpy/isnan.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  typename E::dtype nanmin(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, nanmin);
} // namespace numpy
PYTHONIC_NS_END

#endif
