#ifndef PYTHONIC_INCLUDE_NUMPY_MEDIAN_HPP
#define PYTHONIC_INCLUDE_NUMPY_MEDIAN_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  decltype(std::declval<T>() + 1.) median(types::ndarray<T, pS> const &arr, types::none_type = {});

  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value != 1,
                   types::ndarray<decltype(std::declval<T>() + 1.),
                                  types::array_tuple<long, std::tuple_size<pS>::value - 1>>>
  median(types::ndarray<T, pS> const &arr, long axis);

  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 1, decltype(std::declval<T>() + 1.)>
  median(types::ndarray<T, pS> const &arr, long axis);

  NUMPY_EXPR_TO_NDARRAY0_DECL(median);

  DEFINE_FUNCTOR(pythonic::numpy, median);
} // namespace numpy
PYTHONIC_NS_END

#endif
