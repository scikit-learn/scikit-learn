#ifndef PYTHONIC_INCLUDE_NUMPY_EXPAND_DIMS_HPP
#define PYTHONIC_INCLUDE_NUMPY_EXPAND_DIMS_HPP

#include <pythonic/include/types/ndarray.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <typename T>
  types::ndarray<typename T::dtype, types::array_tuple<long, T::value + 1>>
  expand_dims(T const &input, int axis);

  DEFINE_FUNCTOR(pythonic::numpy, expand_dims);
} // namespace numpy
PYTHONIC_NS_END

#endif
