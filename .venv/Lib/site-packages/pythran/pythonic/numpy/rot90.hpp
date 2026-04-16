#ifndef PYTHONIC_NUMPY_ROT90_HPP
#define PYTHONIC_NUMPY_ROT90_HPP

#include "pythonic/include/numpy/rot90.hpp"

#include "pythonic/numpy/copy.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
  rot90(types::ndarray<T, pS> const &expr, int k)
  {
    auto constexpr N = std::tuple_size<pS>::value;
    if (k % 4 == 0)
      return copy(expr);
    types::array_tuple<long, N> shape = sutils::getshape(expr);
    if (k % 4 != 2)
      std::swap(shape[0], shape[1]);
    types::ndarray<T, types::array_tuple<long, N>> out(shape, builtins::None);
    if (k % 4 == 1) {
      for (int i = 0; i < shape[1]; ++i)
        for (int j = 0; j < shape[0]; ++j)
          out[shape[0] - 1 - j][i] = expr[i][j];
    } else if (k % 4 == 2) {
      for (int i = 0; i < shape[1]; ++i)
        for (int j = 0; j < shape[0]; ++j)
          out[shape[0] - 1 - j][shape[1] - 1 - i] = expr[j][i];
    } else {
      for (int i = 0; i < shape[1]; ++i)
        for (int j = 0; j < shape[0]; ++j)
          out[j][shape[1] - 1 - i] = expr[i][j];
    }
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(rot90)
} // namespace numpy
PYTHONIC_NS_END

#endif
