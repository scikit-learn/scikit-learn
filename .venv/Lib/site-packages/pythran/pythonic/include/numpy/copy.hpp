#ifndef PYTHONIC_INCLUDE_NUMPY_COPY_HPP
#define PYTHONIC_INCLUDE_NUMPY_COPY_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  // list case
  template <class E>
  std::enable_if_t<!types::is_array<E>::value && !types::is_dtype<E>::value,
                   types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>>
  copy(E const &v);

  // scalar / complex case
  template <class E>
  auto copy(E const &v) -> std::enable_if_t<types::is_dtype<E>::value, E>;

  // No copy is required for numpy_expr
  template <class E>
  auto copy(E &&v) -> std::enable_if_t<types::is_array<E>::value, decltype(std::forward<E>(v))>;

  // ndarray case
  template <class T, class pS>
  types::ndarray<T, pS> copy(types::ndarray<T, pS> const &a);

  // transposed ndarray case
  template <class T, class pS>
  types::numpy_texpr<types::ndarray<T, pS>>
  copy(types::numpy_texpr<types::ndarray<T, pS>> const &a);

  DEFINE_FUNCTOR(pythonic::numpy, copy);
} // namespace numpy
PYTHONIC_NS_END

#endif
