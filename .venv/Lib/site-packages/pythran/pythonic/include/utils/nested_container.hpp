#ifndef PYTHONIC_INCLUDE_UTILS_NESTED_CONTAINER_HPP
#define PYTHONIC_INCLUDE_UTILS_NESTED_CONTAINER_HPP

#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include <limits>

PYTHONIC_NS_BEGIN
namespace types
{
  template <class T, class S>
  class sliced_list;
  template <class T>
  class list;
  template <class T, class S>
  class sliced_array;
  template <class T>
  class array;
  template <class T, size_t N, class V>
  struct array_base;
  template <class T>
  struct dynamic_tuple;
} // namespace types

namespace utils
{

  /* compute nested container depth && memory size*/
  template <class T, bool IsArray>
  struct nested_container_depth_helper;

  template <class T>
  struct nested_container_depth_helper<T, false> {
    static const int value = 0;
  };

  template <class T>
  struct nested_container_depth_helper<T, true> {
    static const int value = T::value;
  };

  template <class T>
  struct nested_container_depth {
    static const int value = nested_container_depth_helper<T, types::is_array<T>::value>::value;
  };

  template <class T>
  struct nested_container_depth<types::list<T>> {
    static const int value = 1 + nested_container_depth<T>::value;
  };

  template <class T, class S>
  struct nested_container_depth<types::sliced_list<T, S>> {
    static const int value = 1 + nested_container_depth<T>::value;
  };

  template <class T>
  struct nested_container_depth<types::array<T>> {
    static const int value = 1;
  };

  template <class T, class S>
  struct nested_container_depth<types::sliced_array<T, S>> {
    static const int value = 1;
  };

  template <class T>
  struct nested_container_depth<types::dynamic_tuple<T>> {
    static const int value = 1 + nested_container_depth<T>::value;
  };

  template <class T, size_t N, class V>
  struct nested_container_depth<types::array_base<T, N, V>> {
    static const int value = 1 + nested_container_depth<T>::value;
  };

  template <class T, class sP>
  struct nested_container_depth<types::ndarray<T, sP>> {
    static const int value = std::tuple_size<sP>::value;
  };

  /* Get the size of a container, using recursion on inner container if any
   * FIXME: should be a constexpr?
   * FIXME: why a class and not a function?
   */
  template <class T>
  struct nested_container_size {
    using Type = std::remove_cv_t<std::remove_reference_t<T>>;
    static long flat_size(T const &t);
  };

  /* Recursion stops on bool */
  template <>
  struct nested_container_size<bool> {
    template <class F>
    constexpr static long flat_size(F);
  };

  /* Statically define (by recursion) the type of element inside nested
   * containers */
  template <class T, bool IsArray>
  struct nested_container_value_type_helper;

  template <class T>
  struct nested_container_value_type_helper<T, false> {
    using type = T;
  };

  template <class T>
  struct nested_container_value_type_helper<T, true> {
    using type = typename T::dtype;
  };

  template <class T>
  struct nested_container_value_type {
    using type = typename nested_container_value_type_helper<T, types::is_array<T>::value>::type;
  };

  template <class T>
  struct nested_container_value_type<types::dynamic_tuple<T>> {
    using type = typename nested_container_value_type<T>::type;
  };

  template <class T>
  struct nested_container_value_type<types::list<T>> {
    using type = typename nested_container_value_type<T>::type;
  };

  template <class T, class S>
  struct nested_container_value_type<types::sliced_list<T, S>> {
    using type = typename nested_container_value_type<T>::type;
  };

  template <class T>
  struct nested_container_value_type<types::array<T>> {
    using type = T;
  };

  template <class T, class S>
  struct nested_container_value_type<types::sliced_array<T, S>> {
    using type = T;
  };

  template <class T, size_t N, class V>
  struct nested_container_value_type<types::array_base<T, N, V>> {
    using type = typename nested_container_value_type<T>::type;
  };

  template <class T, class sP>
  struct nested_container_value_type<types::ndarray<T, sP>> {
    using type = T;
  };
} // namespace utils
PYTHONIC_NS_END

#endif
