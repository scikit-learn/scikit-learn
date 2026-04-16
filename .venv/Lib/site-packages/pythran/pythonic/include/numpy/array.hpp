#ifndef PYTHONIC_INCLUDE_NUMPY_ARRAY_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARRAY_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/nested_container.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class dtype = types::dtype_t<typename std::decay_t<T>::dtype>>
  std::enable_if_t<
      types::has_size<std::decay_t<T>>::value,
      types::ndarray<typename dtype::type, types::array_tuple<long, std::decay_t<T>::value>>>
  array(T &&iterable, dtype d = dtype());
  template <class T, class dtype = types::dtype_t<typename std::decay_t<T>::dtype>>
  std::enable_if_t<
      !types::has_size<std::decay_t<T>>::value && !types::is_dtype<std::decay_t<T>>::value,
      types::ndarray<typename dtype::type, types::array_tuple<long, std::decay_t<T>::value>>>
  array(T &&iterable, dtype d = dtype());

  template <class T, class dtype = types::dtype_t<typename types::dtype_of<std::decay_t<T>>::type>>
  std::enable_if_t<!types::has_size<std::decay_t<T>>::value &&
                       types::is_dtype<std::decay_t<T>>::value,
                   typename dtype::type>
  array(T &&non_iterable, dtype d = dtype());

  template <class dtype>
  types::ndarray<typename dtype::type, types::pshape<std::integral_constant<long, 0>>>
      array(std::tuple<>, dtype);

  template <class T, class pS>
  types::ndarray<T, pS> array(types::ndarray<T, pS> const &arr);

  template <class T, size_t N, class V, class dtype = types::dtype_of<T>>
  types::ndarray<typename dtype::type, typename types::array_base<T, N, V>::shape_t>
  array(types::array_base<T, N, V> const &, dtype d = dtype());

  template <class T, size_t N, class V, class dtype = types::dtype_of<T>>
  types::ndarray<typename dtype::type, typename types::array_base<T, N, V>::shape_t>
  array(types::array_base<T, N, V> &&, dtype d = dtype());

  template <class... Ts>
  auto array(std::tuple<Ts...> t)
      -> decltype(array(types::to_array<typename __combined<Ts...>::type>(t)))
  {
    return array(types::to_array<typename __combined<Ts...>::type>(t));
  }

  DEFINE_FUNCTOR(pythonic::numpy, array);
} // namespace numpy
PYTHONIC_NS_END

#endif
