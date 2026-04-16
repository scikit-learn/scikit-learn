#ifndef PYTHONIC_INCLUDE_NUMPY_DTYPE_TYPE_HPP
#define PYTHONIC_INCLUDE_NUMPY_DTYPE_TYPE_HPP
#include "pythonic/include/utils/functor.hpp"
PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace dtype
  {
    template <class T, class V>
    auto type(T const &t, V const &v) -> decltype(t(v));
    DEFINE_FUNCTOR(pythonic::numpy::dtype, type);
  } // namespace dtype
} // namespace numpy

PYTHONIC_NS_END

#endif
