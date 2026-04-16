#ifndef PYTHONIC_NUMPY_DTYPE_TYPE_HPP
#define PYTHONIC_NUMPY_DTYPE_TYPE_HPP

#include "pythonic/include/numpy/dtype/type.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace dtype
  {
    template <class T, class V>
    auto type(T const &t, V const &v) -> decltype(t(v))
    {
      return t(v);
    }
  } // namespace dtype
} // namespace numpy

PYTHONIC_NS_END

#endif
