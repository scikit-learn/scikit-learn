#ifndef PYTHONIC_INCLUDE_TYPES_COMBINED_HPP
#define PYTHONIC_INCLUDE_TYPES_COMBINED_HPP

#include <cstdint>

#include "pythonic/include/types/traits.hpp"
PYTHONIC_NS_BEGIN
namespace types
{
  template <class... Types>
  struct variant_functor;
}
PYTHONIC_NS_END

/* type inference stuff
 */

template <class... Types>
struct __combined;

template <class T>
struct __combined<T> {
  using type = T;
};

template <class T0, class T1, class T2, class... Types>
struct __combined<T0, T1, T2, Types...>
    : __combined<typename __combined<T0, T1>::type, T2, Types...> {
};

template <class T0, class T1>
struct __combined<T0, T1> {
  // callable -> functor
  template <class F0, class F1>
  static pythonic::types::variant_functor<F0, F1> get(std::integral_constant<bool, true>);

  // operator+ exists -> deduce type
  template <class F0, class F1>
  static decltype(std::declval<F0>() + std::declval<F1>()) get(std::integral_constant<bool, false>);

  // operator+ does not exists -> pick first one, better than error
  // note that this is needed because broadcasting is too complex to be modeled
  // by our clumsy type inference scheme
  // so we sometime endup with __combined<indexable_container<...>, int> which
  // only makes sense when broadcasting
  // fortunately, broadcasting is only supported by ndarray, and we already
  // ignore __combined for ndarray
  // so the only thing to do in such situations is « not throw an error »
  template <class F0, class F1>
  static F0 get(...);

  using type = std::conditional_t<
      std::is_same<T0, T1>::value, T0,
      decltype(get<T0, T1>(
          std::integral_constant<bool, pythonic::types::is_callable<T0>::value &&
                                           pythonic::types::is_callable<T1>::value>()))>;
};

template <class T0, class T1>
struct __combined<const T0, T1> : std::add_const<typename __combined<T0, T1>::type> {
};

template <class T0, class T1>
struct __combined<T0, const T1> : std::add_const<typename __combined<T0, T1>::type> {
};

template <class T0, class T1>
struct __combined<T0 &, T1> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 &&, T1> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 const &, T1> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0, T1 &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0, T1 &&> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0, T1 const &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<const T0, T1 const &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<const T0, T1 &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<const T0, T1 &&> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 &, T1 const> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 &&, T1 const> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 const &, T1 const> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 &, T1 const &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 &&, T1 const &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 const &, T1 &> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 const &, T1 &&> : __combined<T0, T1> {
};

template <class T0, class T1>
struct __combined<T0 &, T1 &> : std::add_lvalue_reference<typename __combined<T0, T1>::type> {
};

template <class T0, class T1>
struct __combined<T0 &&, T1 &&> : std::add_rvalue_reference<typename __combined<T0, T1>::type> {
};

template <class T0, class T1>
struct __combined<const T0, const T1> : std::add_const<typename __combined<T0, T1>::type> {
};

template <class T0, class T1>
struct __combined<const T0 &, const T1 &>
    : std::add_lvalue_reference<typename std::add_const<typename __combined<T0, T1>::type>::type> {
};

template <class T>
class container
{
public:
  using value_type = std::remove_cv_t<std::remove_reference_t<T>>;

private:
  container();
};

namespace std
{
  template <size_t I, class T>
  struct tuple_element<I, container<T>> {
    using type = typename container<T>::value_type;
  };
} // namespace std

template <class K, class V>
class indexable_container
{
public:
  using key_type = std::remove_cv_t<std::remove_reference_t<K>>;
  using value_type = std::remove_cv_t<std::remove_reference_t<V>>;

private:
  indexable_container();
};

namespace std
{
  template <size_t I, class K, class V>
  struct tuple_element<I, indexable_container<K, V>> {
    using type = typename indexable_container<K, V>::value_type;
  };
} // namespace std

template <class T>
class dict_container
{
public:
  using value_type = std::remove_cv_t<std::remove_reference_t<T>>;

private:
  dict_container();
};

namespace std
{
  template <size_t I, class T>
  struct tuple_element<I, dict_container<T>> {
    using type = typename dict_container<T>::value_type;
  };
} // namespace std

template <class T>
class indexable
{
public:
  using type = std::remove_cv_t<std::remove_reference_t<T>>;

private:
  indexable();
};

template <class T>
class indexable_dict
{
public:
  using type = std::remove_cv_t<std::remove_reference_t<T>>;

private:
  indexable_dict();
};

template <class K0, class V0, class K1, class V1>
struct __combined<indexable_container<K0, V0>, indexable_container<K1, V1>> {
  using type =
      indexable_container<typename __combined<K0, K1>::type, typename __combined<V0, V1>::type>;
};

template <class K, class V>
struct __combined<indexable<K>, indexable<V>> {
  using type = indexable<typename __combined<K, V>::type>;
};

template <class K, class V>
struct __combined<indexable<K>, container<V>> {
  using type = indexable_container<K, V>;
};

template <class V, class K>
struct __combined<container<V>, indexable<K>> {
  using type = indexable_container<K, V>;
};

template <class K, class V, class W>
struct __combined<indexable_container<K, V>, container<W>> {
  using type = indexable_container<K, typename __combined<V, W>::type>;
};

template <class V, class K, class W>
struct __combined<container<W>, indexable_container<K, V>> {
  using type = indexable_container<K, typename __combined<V, W>::type>;
};

template <class K1, class V1, class K2>
struct __combined<indexable_container<K1, V1>, indexable<K2>> {
  using type = indexable_container<typename __combined<K1, K2>::type, V1>;
};

template <class K1, class V1, class K2>
struct __combined<indexable<K2>, indexable_container<K1, V1>> {
  using type = indexable_container<typename __combined<K1, K2>::type, V1>;
};

template <class A, class B>
struct __combined<container<A>, container<B>> {
  using type = container<typename __combined<A, B>::type>;
};
/* special handling for functors
 * as it's based on a trait, template specialization cannot be used
 * so we rely on operator+ specialization
 * { */

template <class T, class... Types>
struct __combined<T, pythonic::types::variant_functor<Types...>> {
  using type = pythonic::types::variant_functor<T, Types...>;
};

template <class T, class... Types>
struct __combined<pythonic::types::variant_functor<Types...>, T> {
  using type = pythonic::types::variant_functor<T, Types...>;
};

template <class... Types0, class... Types1>
struct __combined<pythonic::types::variant_functor<Types0...>,
                  pythonic::types::variant_functor<Types1...>> {
  using type = pythonic::types::variant_functor<Types0..., Types1...>;
};

/* } */

/* mimic numpy behavior { */
#define SCALAR_COMBINER(Type)                                                                      \
  template <>                                                                                      \
  struct __combined<Type, Type> {                                                                  \
    using type = Type;                                                                             \
  };
SCALAR_COMBINER(bool)
SCALAR_COMBINER(uint8_t)
SCALAR_COMBINER(int8_t)
SCALAR_COMBINER(uint16_t)
SCALAR_COMBINER(int16_t)
SCALAR_COMBINER(uint32_t)
SCALAR_COMBINER(int32_t)
SCALAR_COMBINER(uint64_t)
SCALAR_COMBINER(int64_t)
#undef SCALAR_COMBINER

#endif
