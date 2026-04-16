#ifndef PYTHONIC_INCLUDE_NUMPY_NANARGMIN_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANARGMIN_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  long nanargmin(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, nanargmin);
} // namespace numpy
PYTHONIC_NS_END

#endif
