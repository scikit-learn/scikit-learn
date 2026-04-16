#ifndef PYTHONIC_INCLUDE_TYPES_DICT_HPP
#define PYTHONIC_INCLUDE_TYPES_DICT_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/empty_iterator.hpp"
#include "pythonic/include/types/tuple.hpp"

#include "pythonic/include/utils/allocate.hpp"
#include "pythonic/include/utils/iterator.hpp"
#include "pythonic/include/utils/reserve.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include "pythonic/include/builtins/None.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <utility>

PYTHONIC_NS_BEGIN

namespace types
{

  static const size_t DEFAULT_DICT_CAPACITY = 64;

  struct empty_dict;

  template <class I>
  struct item_iterator_adaptator {
    I base;
    using difference_type = typename std::iterator_traits<I>::difference_type;
    using value_type =
        make_tuple_t<std::remove_cv_t<typename std::iterator_traits<I>::value_type::first_type>,
                     typename std::iterator_traits<I>::value_type::second_type>;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = typename std::iterator_traits<I>::iterator_category;
    item_iterator_adaptator() = default;
    item_iterator_adaptator(I const &i);
    value_type operator*() const;
    bool operator==(item_iterator_adaptator const &other) const
    {
      return base == other.base;
    }
    bool operator!=(item_iterator_adaptator const &other) const
    {
      return base != other.base;
    }
    item_iterator_adaptator &operator++()
    {
      ++base;
      return *this;
    }
    auto operator->() -> decltype(&*base)
    {
      return &*base;
    }
  };

  template <class I>
  struct key_iterator_adaptator {
    I base;
    using difference_type = typename std::iterator_traits<I>::difference_type;
    using value_type = typename std::iterator_traits<I>::value_type::first_type;
    using pointer = typename std::iterator_traits<I>::value_type::first_type *;
    using reference = typename std::iterator_traits<I>::value_type::first_type &;
    using iterator_category = typename std::iterator_traits<I>::iterator_category;
    key_iterator_adaptator() = default;
    key_iterator_adaptator(I const &i);
    value_type operator*() const;
    bool operator==(key_iterator_adaptator const &other) const
    {
      return base == other.base;
    }
    bool operator!=(key_iterator_adaptator const &other) const
    {
      return base != other.base;
    }
    key_iterator_adaptator &operator++()
    {
      ++base;
      return *this;
    }
    auto operator->() -> decltype(&*base)
    {
      return &*base;
    }
  };

  template <class I>
  struct value_iterator_adaptator {
    I base;
    using difference_type = typename std::iterator_traits<I>::difference_type;
    using value_type = typename std::iterator_traits<I>::value_type::second_type;
    using pointer = typename std::iterator_traits<I>::value_type::second_type *;
    using reference = typename std::iterator_traits<I>::value_type::second_type &;
    using iterator_category = typename std::iterator_traits<I>::iterator_category;
    value_iterator_adaptator() = default;
    value_iterator_adaptator(I const &i);
    value_type operator*() const;
    bool operator==(value_iterator_adaptator const &other) const
    {
      return base == other.base;
    }
    bool operator!=(value_iterator_adaptator const &other) const
    {
      return base != other.base;
    }
    value_iterator_adaptator &operator++()
    {
      ++base;
      return *this;
    }
    auto operator->() -> decltype(&*base)
    {
      return &*base;
    }
  };

  template <class D>
  struct dict_items {
    using iterator = typename D::item_const_iterator;
    using value_type = typename iterator::value_type;
    D data;
    dict_items();
    dict_items(D const &d);
    iterator begin() const;
    iterator end() const;
    long size() const;
  };

  template <class D>
  struct dict_keys {
    using iterator = typename D::key_const_iterator;
    using value_type = typename iterator::value_type;
    D data;
    dict_keys();
    dict_keys(D const &d);
    iterator begin() const;
    iterator end() const;
    long size() const;
  };

  template <class D>
  struct dict_values {
    using iterator = typename D::value_const_iterator;
    using value_type = typename iterator::value_type;
    D data;
    dict_values();
    dict_values(D const &d);
    iterator begin() const;
    iterator end() const;
    long size() const;
  };

  // specialization of a map whose key are none_type
  template <class V>
  struct none_type_map {
    using key_type = none_type;
    using mapped_type = V;
    using value_type = std::pair<const key_type, mapped_type>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using iterator = value_type *;
    using const_iterator = const value_type *;

    value_type data_[1];
    bool empty_;

    none_type_map(size_type) : data_(), empty_(true)
    {
    }
    template <class B, class E>
    none_type_map(B begin, E end)
        : data_{begin == end ? value_type() : *begin}, empty_(begin == end)
    {
    }

    mapped_type &operator[](none_type)
    {
      empty_ = false;
      return data_[0].second;
    }

    iterator find(none_type)
    {
      return data_ + empty_;
    }

    const_iterator find(none_type) const
    {
      return data_ + empty_;
    }

    iterator begin()
    {
      return data_ + empty_;
    }
    const_iterator begin() const
    {
      return data_ + empty_;
    }
    iterator end()
    {
      return data_ + 1;
    }
    const_iterator end() const
    {
      return data_ + 1;
    }

    bool empty() const
    {
      return empty_;
    }
    void clear()
    {
      empty_ = true;
      data_[0].second = {};
    }
    const_iterator erase(const_iterator pos)
    {
      clear();
      return end();
    }
    bool erase(none_type)
    {
      bool res = empty_;
      clear();
      return !res;
    }

    size_t size() const
    {
      return empty_ ? 0 : 1;
    }
  };

  template <class K, class V>
  class dict
  {

    // data holder
    using _key_type = std::remove_cv_t<std::remove_reference_t<K>>;
    using _value_type = std::remove_cv_t<std::remove_reference_t<V>>;
    using container_type = std::conditional_t<
        std::is_same<K, none_type>::value, none_type_map<_value_type>,
        std::unordered_map<_key_type, _value_type, std::hash<_key_type>, std::equal_to<_key_type>,
                           utils::allocator<std::pair<const _key_type, _value_type>>>>;

    utils::shared_ref<container_type> data;
    template <class Kp, class Vp>
    friend class dict;

  public:
    // types
    using reference = typename container_type::reference;
    using const_reference = typename container_type::const_reference;
    using iterator =
        utils::comparable_iterator<key_iterator_adaptator<typename container_type::iterator>>;
    using const_iterator =
        utils::comparable_iterator<key_iterator_adaptator<typename container_type::const_iterator>>;
    using item_iterator =
        utils::comparable_iterator<item_iterator_adaptator<typename container_type::iterator>>;
    using item_const_iterator = utils::comparable_iterator<
        item_iterator_adaptator<typename container_type::const_iterator>>;
    using key_iterator =
        utils::comparable_iterator<key_iterator_adaptator<typename container_type::iterator>>;
    using key_const_iterator =
        utils::comparable_iterator<key_iterator_adaptator<typename container_type::const_iterator>>;
    using value_iterator =
        utils::comparable_iterator<value_iterator_adaptator<typename container_type::iterator>>;
    using value_const_iterator = utils::comparable_iterator<
        value_iterator_adaptator<typename container_type::const_iterator>>;
    using size_type = typename container_type::size_type;
    using difference_type = typename container_type::difference_type;
    using value_type = typename container_type::value_type;
    using pointer = typename container_type::pointer;
    using const_pointer = typename container_type::const_pointer;

    // constructors
    dict();
    dict(empty_dict const &);
    dict(std::initializer_list<value_type> l);
    dict(dict<K, V> const &other);
    template <class Kp, class Vp>
    dict(dict<Kp, Vp> const &other);
    template <class B, class E>
    dict(B begin, E end);

    // iterators
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    item_iterator item_begin();
    item_const_iterator item_begin() const;
    item_iterator item_end();
    item_const_iterator item_end() const;
    key_iterator key_begin();
    key_const_iterator key_begin() const;
    key_iterator key_end();
    key_const_iterator key_end() const;
    value_iterator value_begin();
    value_const_iterator value_begin() const;
    value_iterator value_end();
    value_const_iterator value_end() const;

    // dict interface
    operator bool() const;
    V &operator[](K const &key) &;

    template <class OtherKey>
    V &operator[](OtherKey const &key) &
    {
      return (*this)[K(key)];
    }
    V &operator[](K const &key) const &;

    template <class OtherKey>
    V &operator[](OtherKey const &key) const &
    {
      return (*this)[K(key)];
    }

    V &fast(K const &key) &;
    V &fast(K const &key) const &;

    item_const_iterator find(K const &key) const;

    void clear();

    dict<K, V> copy() const;

    template <class W>
    typename __combined<V, W>::type get(K const &key, W d) const;

    none<V> get(K const &key) const;

    template <class W>
    V &setdefault(K const &key, W d);

    none<V> &setdefault(K const &key);

    template <class K0, class W0>
    void update(dict<K0, W0> const &d);

    template <class Iterable>
    void update(Iterable const &d);

    template <class W>
    typename __combined<V, W>::type pop(K const &key, W d);

    V pop(K const &key);

    make_tuple_t<K, V> popitem();

    long size() const;

    dict_items<dict<K, V>> items() const;
    dict_keys<dict<K, V>> keys() const;
    dict_values<dict<K, V>> values() const;

    // type inference stuff
    template <class K_, class V_>
    dict<typename __combined<K, K_>::type, typename __combined<V, V_>::type>
    operator+(dict<K_, V_> const &);

    // id interface
    intptr_t id() const;

    template <class T>
    bool contains(T const &key) const;
  };

  template <class K, class V>
  constexpr dict<K, V> const &as_const(dict<K, V> &t) noexcept
  {
    return t;
  }

  struct empty_dict {

    using value_type = void;
    using iterator = empty_iterator;
    using const_iterator = empty_iterator;

    template <class K, class V>
    dict<K, V> operator+(dict<K, V> const &s);

    empty_dict operator+(empty_dict const &);
    operator bool() const;
    iterator begin() const;
    iterator end() const;
    template <class V>
    bool contains(V const &) const;
  };

  template <class K, class V>
  dict<K, V> operator+(dict<K, V> const &d, empty_dict);
} // namespace types

template <class K, class V>
struct assignable<types::dict<K, V>> {
  using type = types::dict<typename assignable<K>::type, typename assignable<V>::type>;
};

std::ostream &operator<<(std::ostream &os, types::empty_dict const &);

template <class K, class V>
std::ostream &operator<<(std::ostream &os, std::pair<K, V> const &p);

template <class K, class V>
std::ostream &operator<<(std::ostream &os, types::dict<K, V> const &v);
PYTHONIC_NS_END

/* overload std::get */
namespace std
{
  template <size_t I, class K, class V>
  auto get(pythonic::types::dict<K, V> &d) -> decltype(d[I]);

  template <size_t I, class K, class V>
  auto get(pythonic::types::dict<K, V> const &d) -> decltype(d[I]);

  template <size_t I, class K, class V>
  struct tuple_element<I, pythonic::types::dict<K, V>> {
    using type = V;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"
#include "pythonic/include/types/list.hpp"

template <class A>
struct __combined<container<A>, pythonic::types::empty_dict> {
  using type = dict_container<A>;
};

template <class A>
struct __combined<pythonic::types::empty_dict, container<A>> {
  using type = dict_container<A>;
};

template <class A, class B, class C>
struct __combined<container<A>, pythonic::types::dict<C, B>> {
  using type = pythonic::types::dict<C, typename __combined<A, B>::type>;
};

template <class A, class B, class C>
struct __combined<pythonic::types::dict<C, B>, container<A>> {
  using type = pythonic::types::dict<C, typename __combined<A, B>::type>;
};

template <class T>
struct __combined<pythonic::types::empty_dict, pythonic::types::list<T>> {
  using type = pythonic::types::dict<std::tuple_element_t<0, T>, std::tuple_element_t<1, T>>;
};

template <class T, size_t N>
struct __combined<pythonic::types::empty_dict, pythonic::types::static_list<T, N>> {
  using type = pythonic::types::dict<std::tuple_element_t<0, T>, std::tuple_element_t<1, T>>;
};

template <class T>
struct __combined<pythonic::types::list<T>, pythonic::types::empty_dict> {
  using type = pythonic::types::dict<std::tuple_element_t<0, T>, std::tuple_element_t<1, T>>;
};
template <class T, size_t N>
struct __combined<pythonic::types::static_list<T, N>, pythonic::types::empty_dict> {
  using type = pythonic::types::dict<std::tuple_element_t<0, T>, std::tuple_element_t<1, T>>;
};

template <class K0, class V0, class T>
struct __combined<pythonic::types::dict<K0, V0>, pythonic::types::list<T>> {
  using type = pythonic::types::dict<typename __combined<K0, std::tuple_element_t<0, T>>::type,
                                     typename __combined<V0, std::tuple_element_t<1, T>>::type>;
};

template <class K0, class V0, class T, size_t N>
struct __combined<pythonic::types::dict<K0, V0>, pythonic::types::static_list<T, N>> {
  using type = pythonic::types::dict<typename __combined<K0, std::tuple_element_t<0, T>>::type,
                                     typename __combined<V0, std::tuple_element_t<1, T>>::type>;
};

template <class K0, class V0, class T>
struct __combined<pythonic::types::list<T>, pythonic::types::dict<K0, V0>> {
  using type = pythonic::types::dict<typename __combined<K0, std::tuple_element_t<0, T>>::type,
                                     typename __combined<V0, std::tuple_element_t<1, T>>::type>;
};

template <class K0, class V0, class T, size_t N>
struct __combined<pythonic::types::static_list<T, N>, pythonic::types::dict<K0, V0>> {
  using type = pythonic::types::dict<typename __combined<K0, std::tuple_element_t<0, T>>::type,
                                     typename __combined<V0, std::tuple_element_t<1, T>>::type>;
};

template <class K>
struct __combined<indexable<K>, pythonic::types::empty_dict> {
  using type = indexable_dict<K>;
};

template <class K>
struct __combined<pythonic::types::empty_dict, indexable_dict<K>> {
  using type = indexable_dict<K>;
};

template <class K>
struct __combined<indexable_dict<K>, pythonic::types::empty_dict> {
  using type = indexable_dict<K>;
};

template <class K0, class K1, class V1>
struct __combined<pythonic::types::dict<K1, V1>, indexable_dict<K0>> {
  using type = pythonic::types::dict<typename __combined<K0, K1>::type, V1>;
};

template <class K0, class K1, class V1>
struct __combined<indexable_dict<K0>, pythonic::types::dict<K1, V1>> {
  using type = pythonic::types::dict<typename __combined<K0, K1>::type, V1>;
};

template <class K>
struct __combined<pythonic::types::empty_dict, indexable<K>> {
  using type = indexable_dict<K>;
};

template <class K0, class V, class K1>
struct __combined<pythonic::types::dict<K0, V>, indexable<K1>> {
  using type = pythonic::types::dict<typename __combined<K0, K1>::type, V>;
};

template <class K0, class V, class K1>
struct __combined<indexable<K1>, pythonic::types::dict<K0, V>> {
  using type = pythonic::types::dict<typename __combined<K0, K1>::type, V>;
};

template <class K, class V>
struct __combined<pythonic::types::empty_dict, indexable_container<K, V>> {
  using type = pythonic::types::dict<K, V>;
};

template <class K0, class V0, class K1, class V1>
struct __combined<pythonic::types::dict<K0, V0>, indexable_container<K1, V1>> {
  using type =
      pythonic::types::dict<typename __combined<K0, K1>::type, typename __combined<V0, V1>::type>;
};

template <class K0, class V0, class K1, class V1>
struct __combined<indexable_container<K1, V1>, pythonic::types::dict<K0, V0>> {
  using type =
      pythonic::types::dict<typename __combined<K0, K1>::type, typename __combined<V0, V1>::type>;
};

template <class K, class V>
struct __combined<indexable_container<K, V>, pythonic::types::empty_dict> {
  using type = pythonic::types::dict<K, V>;
};

template <class K, class V>
struct __combined<indexable<K>, dict_container<V>> {
  using type = pythonic::types::dict<K, V>;
};

template <class V, class K>
struct __combined<dict_container<V>, indexable<K>> {
  using type = pythonic::types::dict<K, V>;
};

template <class V, class K, class W>
struct __combined<dict_container<V>, indexable_container<K, W>> {
  using type = pythonic::types::dict<K, typename __combined<V, W>::type>;
};

template <class V, class K, class W>
struct __combined<indexable_container<K, W>, dict_container<V>> {
  using type = pythonic::types::dict<K, typename __combined<V, W>::type>;
};

template <class V, class K, class W>
struct __combined<indexable_dict<V>, indexable_container<K, W>> {
  using type = pythonic::types::dict<typename __combined<K, V>::type, W>;
};

template <class V, class K, class W>
struct __combined<indexable_container<K, W>, indexable_dict<V>> {
  using type = pythonic::types::dict<typename __combined<K, V>::type, W>;
};

template <class K, class V, class W>
struct __combined<pythonic::types::dict<K, V>, dict_container<W>> {
  using type = pythonic::types::dict<K, typename __combined<V, W>::type>;
};

template <class V, class K, class W>
struct __combined<dict_container<W>, pythonic::types::dict<K, V>> {
  using type = pythonic::types::dict<K, typename __combined<V, W>::type>;
};

template <class K, class V>
struct __combined<indexable_dict<K>, container<V>> {
  using type = pythonic::types::dict<K, V>;
};

template <class K0, class K1>
struct __combined<indexable_dict<K0>, indexable<K1>> {
  using type = indexable_dict<typename __combined<K0, K1>::type>;
};

template <class K0, class K1>
struct __combined<indexable<K0>, indexable_dict<K1>> {
  using type = indexable_dict<typename __combined<K0, K1>::type>;
};

template <class V, class K>
struct __combined<container<V>, indexable_dict<K>> {
  using type = pythonic::types::dict<K, V>;
};

/* } */
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <typename K, typename V>
struct to_python<types::dict<K, V>> {
  static PyObject *convert(types::dict<K, V> const &v);
};

template <>
struct to_python<types::empty_dict> {
  static PyObject *convert(types::empty_dict);
};

template <typename K, typename V>
struct from_python<types::dict<K, V>> {

  static bool is_convertible(PyObject *obj);

  static types::dict<K, V> convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif

#endif
