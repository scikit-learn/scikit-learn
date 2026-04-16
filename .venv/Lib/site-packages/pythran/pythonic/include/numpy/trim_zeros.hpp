#ifndef PYTHONIC_INCLUDE_NUMPY_TRIMZEROS_HPP
#define PYTHONIC_INCLUDE_NUMPY_TRIMZEROS_HPP

#include "pythonic/include/types/numpy_gexpr.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>

  types::numpy_gexpr<T, types::cstride_normalized_slice<1>>
  trim_zeros(T const &expr, types::str const &trim = "fb");

  DEFINE_FUNCTOR(pythonic::numpy, trim_zeros)
} // namespace numpy
PYTHONIC_NS_END

#endif
