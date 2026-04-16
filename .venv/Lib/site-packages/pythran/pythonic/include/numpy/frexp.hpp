#ifndef PYTHONIC_INCLUDE_NUMPY_FREXP_HPP
#define PYTHONIC_INCLUDE_NUMPY_FREXP_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  std::enable_if_t<std::is_scalar<T>::value, std::tuple<T, int>> frexp(T val);

  template <class E>
  std::enable_if_t<!types::is_dtype<E>::value,
                   std::tuple<types::ndarray<typename E::dtype, typename E::shape_t>,
                              types::ndarray<int, typename E::shape_t>>>
  frexp(E const &arr);

  DEFINE_FUNCTOR(pythonic::numpy, frexp);
} // namespace numpy
PYTHONIC_NS_END

#endif
