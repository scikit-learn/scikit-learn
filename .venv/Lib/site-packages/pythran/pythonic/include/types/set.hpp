#ifndef PYTHONIC_INCLUDE_TYPES_SET_HPP
#define PYTHONIC_INCLUDE_TYPES_SET_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/empty_iterator.hpp"
#include "pythonic/include/types/list.hpp"

#include "pythonic/include/utils/allocate.hpp"
#include "pythonic/include/utils/iterator.hpp"
#include "pythonic/include/utils/reserve.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include "pythonic/include/builtins/in.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <unordered_set>
#include <utility>

PYTHONIC_NS_BEGIN
namespace types
{

  struct empty_set;

  template <class T>
  class set;
} // namespace types
PYTHONIC_NS_END

namespace std
{
  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::set<T>> {
    typedef typename pythonic::types::set<T>::value_type type;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"
template <class A, class B>
struct __combined<container<A>, pythonic::types::set<B>> {
  using type = pythonic::types::set<typename __combined<A, B>::type>;
};

template <class B>
struct __combined<pythonic::types::empty_set, pythonic::types::set<B>> {
  using type = pythonic::types::set<B>;
};

template <class B>
struct __combined<pythonic::types::set<B>, pythonic::types::empty_set> {
  using type = pythonic::types::set<B>;
};

template <class A, class B>
struct __combined<pythonic::types::set<B>, container<A>> {
  using type = pythonic::types::set<typename __combined<A, B>::type>;
};

template <class A, class B>
struct __combined<pythonic::types::list<A>, pythonic::types::set<B>> {
  using type = pythonic::types::set<typename __combined<A, B>::type>;
};

template <class A, class B>
struct __combined<pythonic::types::set<B>, pythonic::types::list<A>> {
  using type = pythonic::types::set<typename __combined<A, B>::type>;
};

template <class A>
struct __combined<pythonic::types::list<A>, pythonic::types::empty_set> {
  using type = pythonic::types::set<A>;
};

template <class A>
struct __combined<pythonic::types::empty_set, pythonic::types::list<A>> {
  using type = pythonic::types::set<A>;
};

template <class K>
struct __combined<indexable<K>, pythonic::types::empty_set> {
  using type = indexable<K>;
};

template <class K>
struct __combined<pythonic::types::empty_set, indexable<K>> {
  using type = indexable<K>;
};

template <class A>
struct __combined<pythonic::types::empty_set, container<A>> {
  using type = pythonic::types::set<A>;
};

template <class K, class V>
struct __combined<indexable<K>, pythonic::types::set<V>> {
  using type = pythonic::types::set<V>;
};

template <class K, class V>
struct __combined<pythonic::types::set<V>, indexable<K>> {
  using type = pythonic::types::set<V>;
};

template <class K, class V1, class V2>
struct __combined<indexable_container<K, V1>, pythonic::types::set<V2>> {
  using type = pythonic::types::set<decltype(std::declval<V1>() + std::declval<V2>())>;
};

template <class K, class V1, class V2>
struct __combined<pythonic::types::set<V2>, indexable_container<K, V1>> {
  using type = pythonic::types::set<decltype(std::declval<V1>() + std::declval<V2>())>;
};

template <class T0, class T1>
struct __combined<pythonic::types::set<T0>, pythonic::types::set<T1>> {
  using type = pythonic::types::set<typename __combined<T0, T1>::type>;
};

/* } */

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T>
  class set
  {

    // data holder
    using _type = std::remove_cv_t<std::remove_reference_t<T>>;
    using container_type =
        std::unordered_set<_type, std::hash<_type>, std::equal_to<_type>, utils::allocator<_type>>;
    utils::shared_ref<container_type> data;

  public:
    template <class U>
    friend class set;

    // types
    using reference = typename container_type::reference;
    using const_reference = typename container_type::const_reference;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using size_type = typename container_type::size_type;
    using difference_type = typename container_type::difference_type;
    using value_type = typename container_type::value_type;
    using allocator_type = typename container_type::allocator_type;
    using pointer = typename container_type::pointer;
    using const_pointer = typename container_type::const_pointer;

    // constructors
    set();
    template <class InputIterator>
    set(InputIterator start, InputIterator stop);
    set(empty_set const &);
    set(std::initializer_list<value_type> l);
    set(set<T> const &other);
    template <class F>
    set(set<F> const &other);

    // iterators
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;

    // modifiers
    T pop();
    void add(const T &x);
    void push_back(const T &x);
    void clear();

    template <class U>
    void discard(U const &elem);

    template <class U>
    void remove(U const &elem);

    // set interface
    operator bool() const;

    long size() const;

    // Misc

    set<T> copy() const;

    template <class U>
    bool isdisjoint(U const &other) const;

    template <class U>
    bool issubset(U const &other) const;

    template <class U>
    bool issuperset(U const &other) const;

    set<T> union_() const;

    template <typename U, typename... Types>
    typename __combined<set<T>, U, Types...>::type union_(U &&other, Types &&...others) const;

    template <typename... Types>
    none_type update(Types &&...others);

    set<T> intersection() const;

    template <typename U, typename... Types>
    typename __combined<set<T>, U, Types...>::type intersection(U const &other,
                                                                Types const &...others) const;

    template <typename... Types>
    void intersection_update(Types const &...others);

    set<T> difference() const;

    template <typename U, typename... Types>
    set<T> difference(U const &other, Types const &...others) const;

    template <class V>
    bool contains(V const &v) const;

    template <typename... Types>
    void difference_update(Types const &...others);

    template <typename U>
    set<typename __combined<T, U>::type> symmetric_difference(set<U> const &other) const;

    template <typename U>
    typename __combined<U, set<T>>::type symmetric_difference(U const &other) const;

    template <typename U>
    void symmetric_difference_update(U const &other);

    // Operators
    template <class U>
    bool operator==(set<U> const &other) const;

    template <class U>
    bool operator<=(set<U> const &other) const;

    template <class U>
    bool operator<(set<U> const &other) const;

    template <class U>
    bool operator>=(set<U> const &other) const;

    template <class U>
    bool operator>(set<U> const &other) const;

    template <class U>
    set<typename __combined<T, U>::type> operator|(set<U> const &other) const;

    template <class U>
    void operator|=(set<U> const &other);

    template <class U>
    set<typename __combined<U, T>::type> operator&(set<U> const &other) const;

    template <class U>
    void operator&=(set<U> const &other);

    template <class U>
    set<T> operator-(set<U> const &other) const;

    template <class U>
    void operator-=(set<U> const &other);

    template <class U>
    set<typename __combined<U, T>::type> operator^(set<U> const &other) const;

    template <class U>
    void operator^=(set<U> const &other);

    intptr_t id() const;

    template <class U>
    friend std::ostream &operator<<(std::ostream &os, set<U> const &v);
  };

  struct empty_set {

    using value_type = void;
    using iterator = empty_iterator;
    using const_iterator = empty_iterator;

    empty_set operator|(empty_set const &);
    template <class T>
    set<T> operator|(set<T> const &s);
    template <class U>
    U operator&(U const &s);
    template <class U>
    U operator-(U const &s);
    empty_set operator^(empty_set const &);
    template <class T>
    set<T> operator^(set<T> const &s);

    template <class... Types>
    none_type update(Types &&...);

    operator bool();
    iterator begin() const;
    iterator end() const;
    template <class V>
    bool contains(V const &) const;

    constexpr long size() const
    {
      return 0;
    }
  };
} // namespace types

template <class T>
struct assignable<types::set<T>> {
  using type = types::set<typename assignable<T>::type>;
};
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <typename T>
struct to_python<types::set<T>> {
  static PyObject *convert(types::set<T> const &v);
};

template <>
struct to_python<types::empty_set> {
  static PyObject *convert(types::empty_set);
};

template <class T>
struct from_python<types::set<T>> {
  static bool is_convertible(PyObject *obj);
  static types::set<T> convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif

#endif
