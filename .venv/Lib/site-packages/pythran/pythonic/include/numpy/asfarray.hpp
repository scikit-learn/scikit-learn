#ifndef PYTHONIC_INCLUDE_NUMPY_ASFARRAY_HPP
#define PYTHONIC_INCLUDE_NUMPY_ASFARRAY_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/numpy/float64.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class dtype = functor::float64>
  auto asfarray(E &&e, dtype d = dtype()) -> decltype(asarray(std::forward<E>(e), d));
  DEFINE_FUNCTOR(pythonic::numpy, asfarray);
} // namespace numpy
PYTHONIC_NS_END

#endif
