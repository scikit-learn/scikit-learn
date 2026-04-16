#ifndef PYTHONIC_INCLUDE_NUMPY_NANARGMAX_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANARGMAX_HPP

#include "pythonic/include/numpy/isnan.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  long nanargmax(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, nanargmax);
} // namespace numpy
PYTHONIC_NS_END

#endif
