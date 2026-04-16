#ifndef PYTHONIC_INCLUDE_NUMPY_TRACE_HPP
#define PYTHONIC_INCLUDE_NUMPY_TRACE_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T>
  typename T::dtype trace(T const &expr, int offset = 0);

  DEFINE_FUNCTOR(pythonic::numpy, trace)
} // namespace numpy
PYTHONIC_NS_END

#endif
