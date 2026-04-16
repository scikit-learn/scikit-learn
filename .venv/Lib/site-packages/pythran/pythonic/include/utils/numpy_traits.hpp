#ifndef PYTHONIC_INCLUDE_UTILS_NUMPY_TRAITS_HPP
#define PYTHONIC_INCLUDE_UTILS_NUMPY_TRAITS_HPP

#include "pythonic/include/types/traits.hpp"

PYTHONIC_NS_BEGIN
namespace types
{

  template <class T, class pS>
  struct ndarray;

  template <class A>
  struct numpy_iexpr;

  template <class A, class F>
  struct numpy_vexpr;

  template <class A, class... S>
  struct numpy_gexpr;

  template <class A>
  struct numpy_texpr;

  template <class A>
  struct numpy_texpr_2;

  template <class O, class... Args>
  struct numpy_expr;

  template <class T>
  class list;

  template <class T, class S>
  class sliced_list;

  struct empty_list;

  template <class T>
  class array;

  template <class T, class S>
  class sliced_array;

  template <class T, size_t N, class V>
  struct array_base;

  template <class T>
  struct dynamic_tuple;

  template <class T>
  struct broadcasted;

  template <class T, class B>
  struct broadcast;

  template <class T>
  struct is_ndarray {
    static constexpr bool value = false;
  };

  template <class T, class pS>
  struct is_ndarray<ndarray<T, pS>> {
    static constexpr bool value = true;
  };

  /* Type trait that checks if a type is a potential numpy expression
   *parameter
   *
   * Only used to write concise expression templates
   */
  template <class T>
  struct is_array {
    static constexpr bool value = false;
  };

  template <class T, class pS>
  struct is_array<ndarray<T, pS>> {
    static constexpr bool value = true;
  };

  template <class A>
  struct is_array<numpy_iexpr<A>> {
    static constexpr bool value = true;
  };

  template <class A, class F>
  struct is_array<numpy_vexpr<A, F>> {
    static constexpr bool value = true;
  };

  template <class A, class... S>
  struct is_array<numpy_gexpr<A, S...>> {
    static constexpr bool value = true;
  };

  template <class A>
  struct is_array<numpy_texpr<A>> {
    static constexpr bool value = true;
  };

  template <class A>
  struct is_array<numpy_texpr_2<A>> {
    static constexpr bool value = true;
  };

  template <class O, class... Args>
  struct is_array<numpy_expr<O, Args...>> {
    static constexpr bool value = true;
  };

  template <class T>
  struct is_numexpr_arg : is_array<T> {
  };

  template <class T>
  struct is_numexpr_arg<T &> : is_numexpr_arg<T> {
  };

  template <class T>
  struct is_numexpr_arg<T const> : is_numexpr_arg<T> {
  };

  template <class T>
  struct is_numexpr_arg<T const &> : is_numexpr_arg<T> {
  };

  template <class T>
  struct is_numexpr_arg<T &&> : is_numexpr_arg<T> {
  };

  template <class T>
  struct is_numexpr_arg<list<T>> {
    static constexpr bool value = is_numexpr_arg<T>::value || is_dtype<T>::value;
  };

  template <class T, class S>
  struct is_numexpr_arg<sliced_list<T, S>> : is_numexpr_arg<list<T>> {
  };

  template <>
  struct is_numexpr_arg<empty_list> {
    static constexpr bool value = true;
  };

  template <class T>
  struct is_numexpr_arg<array<T>> : is_numexpr_arg<list<T>> {
  };

  template <class T, class S>
  struct is_numexpr_arg<sliced_array<T, S>> : is_numexpr_arg<sliced_list<T, S>> {
  };

  template <class T>
  struct is_numexpr_arg<broadcasted<T>> {
    static constexpr bool value = is_numexpr_arg<T>::value || is_dtype<T>::value;
  };

  template <class T, class Tp>
  struct is_numexpr_arg<broadcast<T, Tp>> {
    static constexpr bool value = is_numexpr_arg<T>::value || is_dtype<T>::value;
  };

  template <class T, size_t N, class V>
  struct is_numexpr_arg<array_base<T, N, V>> {
    static constexpr bool value = is_numexpr_arg<T>::value || is_dtype<T>::value;
  };

  template <class T>
  struct is_numexpr_arg<dynamic_tuple<T>> {
    static constexpr bool value = is_numexpr_arg<T>::value || is_dtype<T>::value;
  };

  template <class E>
  struct dtype_of {
    template <class T>
    static typename T::dtype get(typename T::dtype *);
    template <class T>
    static T get(...);
    using type = decltype(get<E>(nullptr));
  };

  template <class T>
  struct has_buffer {
    static constexpr bool value = false;
  };

  template <class T, class pS>
  struct has_buffer<ndarray<T, pS>> {
    static constexpr bool value = true;
  };

  template <class A>
  struct has_buffer<numpy_iexpr<A>> : has_buffer<A> {
  };

  template <class A, class... S>
  struct has_buffer<numpy_gexpr<A, S...>> : has_buffer<A> {
  };

} // namespace types
PYTHONIC_NS_END

#endif
