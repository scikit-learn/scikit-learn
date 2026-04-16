#ifndef PYTHONIC_NUMPY_NDARRAY_ASTYPE_HPP
#define PYTHONIC_NUMPY_NDARRAY_ASTYPE_HPP

#include "pythonic/include/numpy/ndarray/astype.hpp"
#include "pythonic/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {

    template <class E, class dtype>
    auto astype(E &&e, dtype d) -> decltype(asarray(std::forward<E>(e), d))
    {
      return asarray(std::forward<E>(e), d);
    }
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
