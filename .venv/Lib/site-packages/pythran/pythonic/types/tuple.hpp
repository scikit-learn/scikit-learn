#ifndef PYTHONIC_TYPES_TUPLE_HPP
#define PYTHONIC_TYPES_TUPLE_HPP

#include "pythonic/include/types/tuple.hpp"

#include "pythonic/types/assignable.hpp"
#include "pythonic/types/dynamic_tuple.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/nditerator.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/int_.hpp"
#include "pythonic/utils/nested_container.hpp"
#include "pythonic/utils/seq.hpp"

#include <algorithm>
#include <tuple>

namespace std
{
  template <class F0, class S0, class F1, class S1>
  bool operator==(pair<F0, S0> const &self, tuple<F1, S1> const &other)
  {
    return self.first == get<0>(other) && self.second == get<1>(other);
  }
  template <class F0, class S0, class F1, class S1>
  bool operator==(pair<const F0, S0> const &self, tuple<F1, S1> const &other)
  {
    return self.first == get<0>(other) && self.second == get<1>(other);
  }
} // namespace std

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> const &t0,
                                           std::tuple<Types1...> const &t1)
{
  return std::tuple_cat(t0, t1);
}

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> &&t0,
                                           std::tuple<Types1...> const &t1)
{
  return std::tuple_cat(std::forward<Types0...>(t0), t1);
}

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> const &t0,
                                           std::tuple<Types1...> &&t1)
{
  return std::tuple_cat(t0, std::forward<Types1...>(t1));
}

template <class... Types0, class... Types1>
std::tuple<Types0..., Types1...> operator+(std::tuple<Types0...> &&t0, std::tuple<Types1...> &&t1)
{
  return std::tuple_cat(std::forward<Types0...>(t0), std::forward<Types1...>(t1));
}

PYTHONIC_NS_BEGIN

namespace types
{

  /* helper to extract the tail of a tuple, && pop the head
   */

  template <class S, class... Stail>
  std::tuple<Stail...> tuple_tail(std::tuple<S, Stail...> const &t)
  {
    return make_tuple_tail<0>(t, std::make_index_sequence<sizeof...(Stail)>{});
  }

  template <class T, size_t N, class V, class A, size_t... I>
  array_base<T, N, V> array_to_array(A const &a, std::index_sequence<I...>)
  {
    return {(T)std::get<I>(a)...};
  }

  /* inspired by std::array implementation */
  template <typename T, size_t N, class V>
  template <class E>
  long array_base<T, N, V>::_flat_size(E const &e, utils::int_<1>) const
  {
    return N;
  }

  template <typename T, size_t N, class V>
  template <class E, size_t L>
  long array_base<T, N, V>::_flat_size(E const &e, utils::int_<L>) const
  {
    return N * _flat_size(e[0], utils::int_<L - 1>{});
  }

  template <typename T, size_t N, class V>
  long array_base<T, N, V>::flat_size() const
  {
    return _flat_size(*this, utils::int_<value>{});
  }

  template <typename T, size_t N, class V>
  intptr_t array_base<T, N, V>::id() const
  {
    return reinterpret_cast<intptr_t>(&(buffer[0]));
  }

  template <typename T, size_t N, class V>
  void array_base<T, N, V>::fill(const value_type &__u)
  {
    std::fill_n(begin(), size(), __u);
  }

  // Iterators.
  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::iterator array_base<T, N, V>::begin() noexcept
  {
    return {data()};
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_iterator array_base<T, N, V>::begin() const noexcept
  {
    return {data()};
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::iterator array_base<T, N, V>::end() noexcept
  {
    return {data() + N};
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_iterator array_base<T, N, V>::end() const noexcept
  {
    return {data() + N};
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::reverse_iterator array_base<T, N, V>::rbegin() noexcept
  {
    return reverse_iterator(end());
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reverse_iterator array_base<T, N, V>::rbegin() const noexcept
  {
    return const_reverse_iterator(end());
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::reverse_iterator array_base<T, N, V>::rend() noexcept
  {
    return reverse_iterator(begin());
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reverse_iterator array_base<T, N, V>::rend() const noexcept
  {
    return const_reverse_iterator(begin());
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_iterator array_base<T, N, V>::cbegin() const noexcept
  {
    return {&(buffer[0])};
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_iterator array_base<T, N, V>::cend() const noexcept
  {
    return {&(buffer[N])};
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reverse_iterator array_base<T, N, V>::crbegin() const noexcept
  {
    return const_reverse_iterator(end());
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reverse_iterator array_base<T, N, V>::crend() const noexcept
  {
    return const_reverse_iterator(begin());
  }

  // Capacity.
  template <typename T, size_t N, class V>
  constexpr typename array_base<T, N, V>::size_type array_base<T, N, V>::size() const noexcept
  {
    return N;
  }

  template <typename T, size_t N, class V>
  constexpr typename array_base<T, N, V>::size_type array_base<T, N, V>::max_size() const noexcept
  {
    return N;
  }

  template <typename T, size_t N, class V>
  constexpr bool array_base<T, N, V>::empty() const noexcept
  {
    return size() == 0;
  }

  // Element access.
  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::reference array_base<T, N, V>::fast(long n)
  {
    assert(n < (long)size());
    return buffer[n];
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reference array_base<T, N, V>::fast(long n) const noexcept
  {
    assert(n < (long)size());
    return buffer[n];
  }

#ifdef USE_XSIMD
  template <typename T, size_t N, class V>
  template <class vectorizer>
  typename array_base<T, N, V>::simd_iterator array_base<T, N, V>::vbegin(vectorizer) const
  {
    return {&buffer[0]};
  }

  template <typename T, size_t N, class V>
  template <class vectorizer>
  typename array_base<T, N, V>::simd_iterator array_base<T, N, V>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {&buffer[long(size() / vector_size * vector_size)]};
  }
#endif

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::reference array_base<T, N, V>::operator[](long __n)
  {
    auto const index = __n < 0 ? (__n + size()) : __n;
    assert(0 <= index && index < size());
    return buffer[index];
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reference
  array_base<T, N, V>::operator[](long __n) const noexcept
  {
    auto const index = __n < 0 ? (__n + size()) : __n;
    assert(0 <= index && index < size());
    return buffer[index];
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::reference array_base<T, N, V>::front()
  {
    return *begin();
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reference array_base<T, N, V>::front() const
  {
    return *begin();
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::reference array_base<T, N, V>::back()
  {
    return N ? *(end() - 1) : *end();
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_reference array_base<T, N, V>::back() const
  {
    return N ? *(end() - 1) : *end();
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::pointer array_base<T, N, V>::data() noexcept
  {
    return &(buffer[0]);
  }

  template <typename T, size_t N, class V>
  typename array_base<T, N, V>::const_pointer array_base<T, N, V>::data() const noexcept
  {
    return &(buffer[0]);
  }

  template <typename T, size_t N, class V>
  template <size_t M>
  bool array_base<T, N, V>::operator==(array_base<T, M, V> const &other) const
  {
    return N == M && std::equal(begin(), end(), other.begin());
  }

  template <typename T, size_t N, class V>
  template <size_t M>
  bool array_base<T, N, V>::operator!=(array_base<T, M, V> const &other) const
  {
    return !(*this == other);
  }

  template <typename T, size_t N, class V>
  template <size_t M>
  bool array_base<T, N, V>::operator<(array_base<T, M, V> const &other) const
  {
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end());
  }

  template <typename T, size_t N, class V>
  template <class Tp, size_t M>
  array_base<typename __combined<T, Tp>::type, N + M, V>
  array_base<T, N, V>::operator+(array_base<Tp, M, V> const &other) const
  {
    array_base<typename __combined<T, Tp>::type, N + M, V> result;
    auto next = std::copy(begin(), end(), result.begin());
    std::copy(other.begin(), other.end(), next);
    return result;
  }
  template <typename T, size_t N, class V>
  template <class... Types>
  array_base<T, N, V>::operator std::tuple<Types...>() const
  {
    return array_to_tuple(*this, std::make_index_sequence<N>{},
                          typename utils::type_sequence<Types...>{});
  }
  template <typename T, size_t N, class V>
  template <typename Tp>
  array_base<T, N, V>::operator array_base<Tp, N, V>() const
  {
    return array_to_array<Tp, N, V>(*this, std::make_index_sequence<N>{});
  }

  template <typename T, size_t N, class V>
  auto array_base<T, N, V>::to_tuple() const
      -> decltype(array_to_tuple(*this, std::make_index_sequence<N>{},
                                 utils::make_repeated_type<T, N>()))
  {
    return array_to_tuple(*this, std::make_index_sequence<N>{}, utils::make_repeated_type<T, N>());
  }

  template <typename T, size_t N, class V>
  template <class W>
  array_base<T, N, W> array_base<T, N, V>::to_array() const
  {
    return reinterpret_cast<array_base<T, N, W> const &>(*this);
  }

  /* array */
  template <typename T, size_t N, class V>
  std::ostream &operator<<(std::ostream &os, types::array_base<T, N, V> const &v)
  {
    os << "(["[std::is_same<V, types::list_version>::value];
    auto iter = v.begin();
    if (iter != v.end()) {
      while (iter + 1 != v.end())
        os << *iter++ << ", ";
      os << *iter;
    }
    return os << ")]"[std::is_same<V, types::list_version>::value];
  }

  template <class T, size_t N, class V, class... Types>
  auto operator+(std::tuple<Types...> const &t, types::array_base<T, N, V> const &lt)
      -> decltype(std::tuple_cat(t, lt.to_tuple()))
  {
    return std::tuple_cat(t, lt.to_tuple());
  }

  template <class T, size_t N, class V, class... Types>
  auto operator+(types::array_base<T, N, V> const &lt, std::tuple<Types...> const &t)
      -> decltype(std::tuple_cat(lt.to_tuple(), t))
  {
    return std::tuple_cat(lt.to_tuple(), t);
  }

  template <class T, size_t N>
  dynamic_tuple<T> array_base_slicer::operator()(array_tuple<T, N> const &b, slice const &s)
  {
    normalized_slice ns = s.normalize(b.size());
    array_tuple<T, N> tmp;
    for (long j = 0; j < ns.size(); ++j)
      tmp[j] = b[ns.lower + j * ns.step];
    return {&tmp[0], &tmp[ns.size()]};
  }

  template <class T, size_t N, long stride>
  dynamic_tuple<T> array_base_slicer::operator()(array_tuple<T, N> const &b,
                                                 cstride_slice<stride> const &s)
  {
    auto ns = s.normalize(b.size());
    if (stride == 1) {
      return {&b[ns.lower], &b[ns.upper]};
    } else {
      array_tuple<T, N> tmp;
      for (long j = 0; j < ns.size(); ++j)
        tmp[j] = b[ns.lower + j * ns.step];
      return {&tmp[0], &tmp[ns.size()]};
    }
  }

  template <class T, size_t N>
  dynamic_tuple<T> array_base_slicer::operator()(array_tuple<T, N> const &b,
                                                 fast_contiguous_slice const &s)
  {
    auto cns = s.normalize(b.size());
    return {&b[cns.lower], &b[cns.upper]};
  }
} // namespace types
PYTHONIC_NS_END

/* hashable tuples, as proposed in
 * http://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
 */
namespace
{

  inline size_t hash_combiner(size_t left, size_t right) // replacable
  {
    return left ^ right;
  }

  template <size_t index, class... types>
  size_t hash_impl<index, types...>::operator()(size_t a, const std::tuple<types...> &t) const
  {
    using nexttype = std::tuple_element_t<index, std::tuple<types...>>;
    hash_impl<index - 1, types...> next;
    size_t b = std::hash<nexttype>()(std::get<index>(t));
    return next(hash_combiner(a, b), t);
  }

  template <class... types>
  size_t hash_impl<0, types...>::operator()(size_t a, const std::tuple<types...> &t) const
  {
    using nexttype = std::tuple_element_t<0, std::tuple<types...>>;
    size_t b = std::hash<nexttype>()(std::get<0>(t));
    return hash_combiner(a, b);
  }
} // namespace

/* specialize std::hash */
namespace std
{
  template <class... Types>
  size_t hash<std::tuple<Types...>>::operator()(std::tuple<Types...> const &t) const
  {
    const size_t begin = std::tuple_size<std::tuple<Types...>>::value - 1;
    return hash_impl<begin, Types...>()(1, t); // 1 should be some largervalue
  }

  template <typename T, size_t N, class V>
  size_t hash<pythonic::types::array_base<T, N, V>>::operator()(
      pythonic::types::array_base<T, N, V> const &l) const
  {
    size_t seed = 0;
    hash<T> h;
    for (auto const &iter : l)
      seed ^= h(iter) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
} // namespace std

PYTHONIC_NS_BEGIN

namespace types
{
  template <class Tuple, size_t I>
  void print_tuple(std::ostream &os, Tuple const &t, utils::int_<I>)
  {
    print_tuple(os, t, utils::int_<I - 1>());
    os << ", " << std::get<I>(t);
  }

  template <class Tuple>
  void print_tuple(std::ostream &os, Tuple const &t, utils::int_<0>)
  {
    os << std::get<0>(t);
  }
} // namespace types
PYTHONIC_NS_END

namespace std
{
  template <class... Args>
  ostream &operator<<(ostream &os, tuple<Args...> const &t)
  {
    os << '(';
    pythonic::types::print_tuple(os, t, pythonic::utils::int_<sizeof...(Args) - 1>());
    return os << ')';
  }
} // namespace std
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/include/utils/fwd.hpp"
#include "pythonic/include/utils/seq.hpp"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

#ifdef Py_LIMITED_API
#define PyTuple_SET_ITEM PyTuple_SetItem
#define PyList_SET_ITEM PyList_SetItem
#endif

template <typename K, typename V>
PyObject *to_python<std::pair<K, V>>::convert(std::pair<K, V> const &t)
{
  PyObject *out = PyTuple_New(2);
  PyTuple_SET_ITEM(out, 0, ::to_python(std::get<0>(t)));
  PyTuple_SET_ITEM(out, 1, ::to_python(std::get<1>(t)));
  return out;
}

template <typename... Tys>
PyObject *to_python<types::pshape<Tys...>>::convert(types::pshape<Tys...> const &t)
{
  return ::to_python(t.array());
}

template <typename... Types>
template <size_t... S>
PyObject *to_python<std::tuple<Types...>>::

    do_convert(std::tuple<Types...> const &t, std::index_sequence<S...>)
{
  PyObject *out = PyTuple_New(sizeof...(Types));
  (void)std::initializer_list<bool>{
      (PyTuple_SET_ITEM(out, S, ::to_python(std::get<S>(t))), true)...};
  return out;
}

template <typename... Types>
PyObject *to_python<std::tuple<Types...>>::convert(std::tuple<Types...> const &t)
{
  return do_convert(t, std::make_index_sequence<sizeof...(Types)>());
}

template <typename T, size_t N>
template <size_t... S>
PyObject *to_python<types::array_tuple<T, N>>::do_convert(types::array_tuple<T, N> const &t,
                                                          std::index_sequence<S...>)
{
  PyObject *out = PyTuple_New(N);
  (void)std::initializer_list<bool>{
      (PyTuple_SET_ITEM(out, S, ::to_python(std::get<S>(t))), true)...};
  return out;
}

template <typename T, size_t N>
template <size_t... S>
PyObject *to_python<types::static_list<T, N>>::do_convert(types::static_list<T, N> const &t,
                                                          std::index_sequence<S...>)
{
  PyObject *out = PyList_New(N);
  (void)std::initializer_list<bool>{
      (PyList_SET_ITEM(out, S, ::to_python(std::get<S>(t))), true)...};
  return out;
}

template <typename T, size_t N>
PyObject *to_python<types::array_tuple<T, N>>::convert(types::array_tuple<T, N> const &t)
{
  return do_convert(t, std::make_index_sequence<N>());
}

template <typename T, size_t N>
PyObject *to_python<types::static_list<T, N>>::convert(types::static_list<T, N> const &t)
{
  return do_convert(t, std::make_index_sequence<N>());
}

template <typename... Types>
template <size_t... S>
bool from_python<std::tuple<Types...>>

    ::do_is_convertible(PyObject *obj, typename std::index_sequence<S...>)
{
  bool checks[] = {::is_convertible<std::tuple_element_t<S, std::tuple<Types...>>>(
#ifdef Py_LIMITED_API
      PyTuple_GetItem(obj, S)
#else
      PyTuple_GET_ITEM(obj, S)
#endif
          )...};
  return std::find(std::begin(checks), std::end(checks), false) == std::end(checks);
}

template <typename... Types>
bool from_python<std::tuple<Types...>>::is_convertible(PyObject *obj)
{
  if (PyTuple_Check(obj)) {
#ifdef Py_LIMITED_API
    auto n = PyTuple_Size(obj);
#else
    auto n = PyTuple_GET_SIZE(obj);
#endif
    if (n == sizeof...(Types)) {
      return do_is_convertible(obj, std::make_index_sequence<sizeof...(Types)>());
    }
  }
  return false;
}

template <typename... Types>
template <size_t... S>
std::tuple<Types...>
from_python<std::tuple<Types...>>::do_convert(PyObject *obj, typename std::index_sequence<S...>)
{
  return std::tuple<Types...>{::from_python<std::tuple_element_t<S, std::tuple<Types...>>>(
#ifdef Py_LIMITED_API
      PyTuple_GetItem(obj, S)
#else
      PyTuple_GET_ITEM(obj, S)
#endif
          )...};
}
template <typename... Types>
std::tuple<Types...> from_python<std::tuple<Types...>>::convert(PyObject *obj)
{
  return do_convert(obj, std::make_index_sequence<sizeof...(Types)>());
}

template <typename T, size_t N>
bool from_python<types::array_tuple<T, N>>::

    is_convertible(PyObject *obj)
{
  if (PyTuple_Check(obj)) {
#ifdef Py_LIMITED_API
    auto n = PyTuple_Size(obj);
#else
    auto n = PyTuple_GET_SIZE(obj);
#endif
    if (n == N) {
      return ::is_convertible<T>(
#ifdef Py_LIMITED_API
          PyTuple_GetItem(obj, 0)
#else
          PyTuple_GET_ITEM(obj, 0)
#endif
      );
    }
  }
  return false;
}

template <typename T, size_t N>
template <size_t... S>
types::array_tuple<T, N>
from_python<types::array_tuple<T, N>>::do_convert(PyObject *obj, typename std::index_sequence<S...>)
{
  return {::from_python<T>(
#ifdef Py_LIMITED_API
      PyTuple_GetItem(obj, S)
#else
      PyTuple_GET_ITEM(obj, S)
#endif
          )...};
}
template <typename T, size_t N>
types::array_tuple<T, N> from_python<types::array_tuple<T, N>>::

    convert(PyObject *obj)
{
  return do_convert(obj, std::make_index_sequence<N>());
}

#ifdef Py_LIMITED_API
#undef PyTuple_SET_ITEM
#undef PyList_SET_ITEM
#endif
PYTHONIC_NS_END
#endif

#endif
