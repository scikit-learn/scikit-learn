#ifndef PYTHONIC_INCLUDE_NUMPY_BINCOUNT_HPP
#define PYTHONIC_INCLUDE_NUMPY_BINCOUNT_HPP

#include "pythonic/include/numpy/max.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 1, types::ndarray<long, types::pshape<long>>>
  bincount(types::ndarray<T, pS> const &expr, types::none_type weights = builtins::None,
           types::none<long> minlength = builtins::None);

  template <class T, class E, class pS>
  std::enable_if_t<
      std::tuple_size<pS>::value == 1,
      types::ndarray<decltype(std::declval<long>() * std::declval<typename E::dtype>()),
                     types::pshape<long>>>
  bincount(types::ndarray<T, pS> const &expr, E const &weights,
           types::none<long> minlength = builtins::None);

  NUMPY_EXPR_TO_NDARRAY0_DECL(bincount);

  DEFINE_FUNCTOR(pythonic::numpy, bincount);
} // namespace numpy
PYTHONIC_NS_END

#endif
