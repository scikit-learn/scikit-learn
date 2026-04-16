#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_OP_HELPER_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_OP_HELPER_HPP

#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/meta.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T0, class... Types>
  struct all_valid_arg {
    static constexpr bool value = all_valid_arg<T0>::value && all_valid_arg<Types...>::value;
  };

  template <class T>
  struct all_valid_arg<T> {
    static constexpr bool value =
        (is_numexpr_arg<T>::value || is_complex<T>::value || std::is_scalar<T>::value);
  };

  template <class T0, class... Types>
  struct any_numexpr_arg {
    static constexpr bool value = any_numexpr_arg<T0>::value || any_numexpr_arg<Types...>::value;
  };

  template <class T>
  struct any_numexpr_arg<T> : is_numexpr_arg<T> {
  };

  template <class... Types>
  struct valid_numexpr_parameters {
    static constexpr bool value =
        any_numexpr_arg<Types...>::value && all_valid_arg<Types...>::value;
  };

  template <>
  struct valid_numexpr_parameters<> : std::false_type {
  };

  template <class T0, class... Types>
  struct any_numop_arg {
    static constexpr bool value = any_numop_arg<T0>::value || any_numop_arg<Types...>::value;
  };

  template <class T>
  struct any_numop_arg<T> : is_numexpr_arg<T> {
  };
  template <>
  struct any_numop_arg<empty_list> : std::false_type {
  };
  template <class T>
  struct any_numop_arg<dynamic_tuple<T>> : std::false_type {
  };
  template <class T>
  struct any_numop_arg<list<T>> : std::false_type {
  };
  template <class T, class S>
  struct any_numop_arg<sliced_list<T, S>> : std::false_type {
  };
  template <class T>
  struct any_numop_arg<array<T>> : std::false_type {
  };
  template <class T, class S>
  struct any_numop_arg<sliced_array<T, S>> : std::false_type {
  };
  template <class T, size_t N, class V>
  struct any_numop_arg<array_base<T, N, V>> : std::false_type {
  };

  template <class... Types>
  struct valid_numop_parameters {
    static constexpr bool value = any_numop_arg<Types...>::value && all_valid_arg<Types...>::value;
  };

  template <>
  struct valid_numop_parameters<> : std::false_type {
  };

  template <class T0, class T1, bool numexprarg = valid_numexpr_parameters<T0, T1>::value,
            bool T0_number = is_dtype<T0>::value, bool T1_number = is_dtype<T1>::value>
  struct the_common_type {
    using type = std::conditional_t < std::decay_t<T0>::value<std::decay_t<T1>::value, T1, T0>;
  };

  template <class T0, class T1>
  struct the_common_type<T0, T1, true, false, true> {
    using type = T0;
  };

  template <class T0, class T1>
  struct the_common_type<T0, T1, true, true, false> {
    using type = T1;
  };

  template <class T0, class T1, bool b0, bool b1>
  struct the_common_type<T0, T1, false, b0, b1> {
    using type = T0; // keep the first one! It's important for the type
                     // adaptation to avoid type promotion
  };

  template <class... Types>
  struct common_type;

  template <class T0>
  struct common_type<T0> {
    using type = T0;
  };

  template <class T0, class T1>
  struct common_type<T0, T1> {
    using type = typename the_common_type<T0, T1>::type;
  };

  template <class T0, class T1, class... Types>
  struct common_type<T0, T1, Types...> {
    using type = typename common_type<typename common_type<T0, T1>::type,
                                      typename common_type<Types...>::type>::type;
  };

  /* An adapted type creates a type that has the same shape as C and the same
   * dtype as T
   * to the exception of broadcasted constants that may take the dtype of C
   * instead
   */
  template <class T, class C, bool same, bool scalar>
  struct adapted_type;

  template <class T, class C, bool scalar>
  struct adapted_type<T, C, true, scalar> {
    using type = T;
  };

  template <class T, class C>
  struct adapted_type<T, C, false, true> {
    using type = broadcast<typename std::decay_t<C>::dtype, std::decay_t<T>>;
  };

  template <class T, size_t N>
  struct broadcasted_n;
  template <class T>
  struct broadcasted_n<T, 1> {
    using type = broadcasted<T>;
  };
  template <class T, size_t N>
  struct broadcasted_n {
    using type = broadcasted<typename broadcasted_n<T, N - 1>::type>;
  };

  constexpr size_t absdiff(size_t x, size_t y)
  {
    return x > y ? x - y : y - x;
  }

  template <class T, class C>
  struct adapted_type<T, C, false, false> {
    using type = typename broadcasted_n<T, absdiff(std::remove_reference_t<T>::value,
                                                   std::remove_reference_t<C>::value)>::type;
  };

  template <class T, class... OtherTypes>
  struct adapt_type {
    using ctype = typename common_type<std::decay_t<T>, OtherTypes...>::type;
    static constexpr bool isdtype = is_dtype<std::decay_t<T>>::value;
    using type =
        typename adapted_type<T, ctype, std::is_same<std::decay_t<T>, ctype>::value, isdtype>::type;
  };
  template <class T, class Tp, class... OtherTypes>
  struct adapt_type<broadcast<T, Tp>, OtherTypes...> {
    using type = broadcast<T, Tp>;
  };

  /* A reshaped type create a type that has the same shape as C && the same
   * dtype as T
   * To the opposite of an adapted type, it does *not* changes constants type
   */
  template <class T, class C, bool same, bool scalar>
  struct reshaped_type;

  template <class T, class C, bool scalar>
  struct reshaped_type<T, C, true, scalar> {
    using type = T;
  };

  template <class T, class C>
  struct reshaped_type<T, C, false, true> {
    using type = broadcast<std::decay_t<T>, std::decay_t<T>>;
  };

  template <class T, class C>
  struct reshaped_type<T, C, false, false> {
    using type = broadcasted<T>;
  };

  template <class T, class... OtherTypes>
  struct reshape_type {
    using ctype = typename common_type<std::decay_t<T>, OtherTypes...>::type;
    static constexpr bool isdtype = is_dtype<std::decay_t<T>>::value;
    using type = typename reshaped_type<T, ctype, std::is_same<std::decay_t<T>, ctype>::value,
                                        isdtype>::type;
  };
  template <class T>
  struct is_array_index : std::false_type {
  };
  template <size_t N>
  struct is_array_index<array_tuple<long, N>> : std::true_type {
  };
} // namespace types
PYTHONIC_NS_END

#endif
