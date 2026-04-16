#ifndef PYTHONIC_INCLUDE_NUMPY_ATLEAST3D_HPP
#define PYTHONIC_INCLUDE_NUMPY_ATLEAST3D_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  std::enable_if_t<types::is_dtype<T>::value,
                   types::ndarray<T, types::pshape<std::integral_constant<long, 1>,
                                                   std::integral_constant<long, 1>,
                                                   std::integral_constant<long, 1>>>>
  atleast_3d(T t);
  template <class T>
  auto atleast_3d(T const &t) -> std::enable_if_t<
      (!types::is_dtype<T>::value) && (T::value == 1),
      types::ndarray<typename T::dtype, types::pshape<std::integral_constant<long, 1>,
                                                      std::tuple_element_t<0, typename T::shape_t>,
                                                      std::integral_constant<long, 1>>>>;

  template <class T>
  auto atleast_3d(T const &t) -> std::enable_if_t<
      (!types::is_dtype<T>::value) && (T::value == 2),
      types::ndarray<typename T::dtype, types::pshape<std::tuple_element_t<0, typename T::shape_t>,
                                                      std::tuple_element_t<1, typename T::shape_t>,
                                                      std::integral_constant<long, 1>>>>;

  template <class T>
  auto atleast_3d(T const &t)
      -> std::enable_if_t<(!types::is_dtype<T>::value) && T::value >= 3, decltype(asarray(t))>;

  DEFINE_FUNCTOR(pythonic::numpy, atleast_3d);
} // namespace numpy
PYTHONIC_NS_END

#endif
