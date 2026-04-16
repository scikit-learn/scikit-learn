#ifndef PYTHONIC_NUMPY_EXPAND_DIMS_HPP
#define PYTHONIC_NUMPY_EXPAND_DIMS_HPP

#include "pythonic/include/numpy/expand_dims.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <typename T>
  types::ndarray<typename T::dtype, types::array_tuple<long, T::value + 1>>
  expand_dims(T const &input, int axis)
  {
    const long N = T::value;
    if (axis == -1)
      axis += N + 1;
    types::array_tuple<long, N + 1> dim_array;
    auto in_shape = sutils::getshape(input);
    long ii, jj;
    for (ii = jj = 0; ii < N + 1; ii++) {
      if (ii == axis) {
        dim_array[ii] = 1;
      } else {
        dim_array[ii] = in_shape[jj++];
      }
    }

    return numpy::functor::asarray{}(input).reshape(dim_array);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
