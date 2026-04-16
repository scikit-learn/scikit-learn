#ifndef PYTHONIC_NUMPY_TRIU_HPP
#define PYTHONIC_NUMPY_TRIU_HPP

#include "pythonic/include/numpy/triu.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, pS> triu(types::ndarray<T, pS> const &expr, int k)
  {
    types::ndarray<T, pS> out(expr._shape, builtins::None);
    for (int i = 0; i < std::get<0>(expr._shape); ++i)
      for (long j = 0; j < std::get<1>(expr._shape); ++j)
        if (j - i >= k)
          out[i][j] = expr[i][j];
        else
          out[i][j] = 0;
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(triu)
} // namespace numpy
PYTHONIC_NS_END

#endif
