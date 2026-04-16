#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_ASTYPE_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_ASTYPE_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {

    template <class E, class dtype>
    auto astype(E &&e, dtype d) -> decltype(asarray(std::forward<E>(e), d));

    DEFINE_FUNCTOR(pythonic::numpy::ndarray, astype);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
