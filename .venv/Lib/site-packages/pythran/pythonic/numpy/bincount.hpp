#ifndef PYTHONIC_NUMPY_BINCOUNT_HPP
#define PYTHONIC_NUMPY_BINCOUNT_HPP

#include "pythonic/include/numpy/bincount.hpp"

#include "pythonic/numpy/max.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 1, types::ndarray<long, types::pshape<long>>>
  bincount(types::ndarray<T, pS> const &expr, types::none_type weights, types::none<long> minlength)
  {
    long length = 0;
    if (minlength)
      length = (long)minlength;
    length = std::max<long>(length, 1 + max(expr));
    types::ndarray<long, types::pshape<long>> out(types::pshape<long>(length), 0L);
    for (auto iter = expr.fbegin(), end = expr.fend(); iter != end; ++iter)
      ++out[*iter];
    return out;
  }

  template <class T, class E, class pS>
  std::enable_if_t<
      std::tuple_size<pS>::value == 1,
      types::ndarray<decltype(std::declval<long>() * std::declval<typename E::dtype>()),
                     types::pshape<long>>>
  bincount(types::ndarray<T, pS> const &expr, E const &weights, types::none<long> minlength)
  {
    long length = 0;
    if (minlength)
      length = (long)minlength;
    length = std::max<long>(length, 1 + max(expr));
    std::enable_if_t<
        std::tuple_size<pS>::value == 1,
        types::ndarray<decltype(std::declval<long>() * std::declval<typename E::dtype>()),
                       types::pshape<long>>>
        out(types::pshape<long>(length), 0L);
    auto iweight = weights.begin();
    for (auto iter = expr.fbegin(), end = expr.fend(); iter != end; ++iter, ++iweight)
      out[*iter] += *iweight;
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(bincount);
} // namespace numpy
PYTHONIC_NS_END

#endif
