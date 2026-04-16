#ifndef PYTHONIC_INCLUDE_NUMPY_ATLEAST1D_HPP
#define PYTHONIC_INCLUDE_NUMPY_ATLEAST1D_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  std::enable_if_t<types::is_dtype<T>::value,
                   types::ndarray<T, types::pshape<std::integral_constant<long, 1>>>>
  atleast_1d(T t);

  template <class T>
  auto atleast_1d(T const &t) -> std::enable_if_t<!types::is_dtype<T>::value, decltype(asarray(t))>;

  DEFINE_FUNCTOR(pythonic::numpy, atleast_1d);
} // namespace numpy
PYTHONIC_NS_END

#endif
