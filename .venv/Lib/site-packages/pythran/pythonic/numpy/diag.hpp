#ifndef PYTHONIC_NUMPY_DIAG_HPP
#define PYTHONIC_NUMPY_DIAG_HPP

#include "pythonic/include/numpy/diag.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 2, types::ndarray<T, types::pshape<long>>>
  diag(types::ndarray<T, pS> const &a, long k)
  {
    auto &&a_shape = a._shape;
    utils::shared_ref<types::raw_array<T>> buffer(
        std::max(std::get<0>(a_shape), std::get<1>(a_shape)));
    types::pshape<long> shape = 0;
    auto iter = buffer->data;
    if (k >= 0)
      for (int i = 0, j = k; i < std::get<0>(a_shape) && j < std::get<1>(a_shape);
           ++i, ++j, ++std::get<0>(shape))
        *iter++ = a[i][j];
    else
      for (int i = -k, j = 0; i < std::get<0>(a_shape) && j < std::get<1>(a_shape);
           ++i, ++j, ++std::get<0>(shape))
        *iter++ = a[i][j];
    return {buffer, shape};
  }

  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 1, types::ndarray<T, types::array_tuple<long, 2>>>
  diag(types::ndarray<T, pS> const &a, long k)
  {
    long n = a.flat_size() + std::abs(k);
    types::ndarray<T, types::array_tuple<long, 2>> out(types::make_tuple(n, n), 0);
    if (k >= 0)
      for (long i = 0, j = k; i < n && j < n; ++i, ++j)
        out[i][j] = a[i];
    else
      for (long i = -k, j = 0; i < n && j < n; ++i, ++j)
        out[i][j] = a[j];
    return out;
  }

  template <class T>
  auto diag(types::list<T> const &a, long k) -> decltype(diag(asarray(a), k))
  {
    return diag(asarray(a), k);
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(diag);
} // namespace numpy
PYTHONIC_NS_END

#endif
