#ifndef PYTHONIC_NUMPY_ASFARRAY_HPP
#define PYTHONIC_NUMPY_ASFARRAY_HPP

#include "pythonic/include/numpy/asfarray.hpp"
#include "pythonic/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class dtype>
  auto asfarray(E &&e, dtype d) -> decltype(asarray(std::forward<E>(e), d))
  {
    static_assert(std::is_floating_point<typename dtype::type>::value,
                  "expected a floating point type");
    return asarray(std::forward<E>(e), d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
