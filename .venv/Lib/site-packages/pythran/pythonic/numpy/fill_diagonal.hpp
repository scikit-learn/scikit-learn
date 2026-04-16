#ifndef PYTHONIC_NUMPY_FILL_DIAGONAL_HPP
#define PYTHONIC_NUMPY_FILL_DIAGONAL_HPP

#include "pythonic/include/numpy/fill_diagonal.hpp"

#include "pythonic/types/NoneType.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::none_type fill_diagonal(E &&expr, typename std::decay_t<E>::dtype fill_value)
  {
    constexpr auto N = std::decay_t<E>::value;
    types::array_tuple<long, N> indices;
    for (long i = 0, n = sutils::min(expr); i < n; ++i) {
      std::fill(indices.begin(), indices.end(), i);
      expr.fast(indices) = fill_value;
    }
    return {};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
