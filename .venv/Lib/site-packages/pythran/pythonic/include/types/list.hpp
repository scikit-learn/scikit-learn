#ifndef PYTHONIC_INCLUDE_TYPES_LIST_HPP
#define PYTHONIC_INCLUDE_TYPES_LIST_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/empty_iterator.hpp"
#include "pythonic/include/types/nditerator.hpp"
#include "pythonic/include/types/slice.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/types/vectorizable_type.hpp"
#include "pythonic/include/utils/allocate.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/utils/reserve.hpp"
#include "pythonic/include/utils/shared_ref.hpp"

#include <algorithm>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  using container = std::vector<T, utils::allocator<T>>;

  /* forward declaration */
  struct empty_list;
  template <class T>
  class list;
  template <class T, class S>
  class sliced_list;
  template <class T, class pS>
  struct ndarray;
  template <class... Tys>
  struct pshape;

  /* list view */
  template <class T, class S = slice>
  class sliced_list
  {

    // data holder
    using _type = std::remove_cv_t<std::remove_reference_t<T>>;
    typedef container<_type> container_type;
    utils::shared_ref<container_type> _data;

    template <class U>
    friend class list;

    typename S::normalized_type slicing;

  public:
    //  types
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef nditerator<sliced_list> iterator;
    typedef const_nditerator<sliced_list> const_iterator;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::value_type value_type;
    typedef typename container_type::allocator_type allocator_type;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
    typedef typename container_type::reverse_iterator reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;

    // minimal ndarray interface
    typedef typename utils::nested_container_value_type<sliced_list>::type dtype;
    static const size_t value = utils::nested_container_depth<sliced_list>::value;
    static_assert(value != 0, "valid shape");
    static const bool is_vectorizable =
        types::is_vectorizable_dtype<dtype>::value && !std::is_same<S, slice>::value;
    static const bool is_flat = std::is_same<slice, S>::value;
    static const bool is_strided = std::is_same<slice, S>::value;

    using shape_t = types::array_tuple<long, value>;
    template <size_t I>
    auto shape() const -> decltype(details::extract_shape(*this, utils::int_<I>{}))
    {
      return details::extract_shape(*this, utils::int_<I>{});
    }

    // constructor
    sliced_list();
    sliced_list(sliced_list<T, S> const &s);
    sliced_list(list<T> const &other, S const &s);
    template <class Sn>
    sliced_list(utils::shared_ref<container_type> const &other, Sn const &s);

    // assignment
    sliced_list &operator=(list<T> const &);
    sliced_list &operator=(sliced_list<T, S> const &);
    list<T> operator+(list<T> const &) const;
    template <size_t N, class V>
    list<T> operator+(array_base<T, N, V> const &) const;
    template <class Tp, class Sp>
    list<typename __combined<T, Tp>::type> operator+(sliced_list<Tp, Sp> const &) const;

    // iterators
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;

    // size
    long size() const;
    explicit operator bool() const;

    // accessors
    const_reference fast(long i) const;
    const_reference operator[](long i) const;
    reference operator[](long i);
    template <class Sp>
    std::enable_if_t<is_slice<Sp>::value,
                     sliced_list<T, decltype(std::declval<S>() * std::declval<Sp>())>>
    operator[](Sp s) const;

    template <class... Indices>
    dtype load(long index0, long index1, Indices... indices) const
    {
      return fast(index0).load(index1, indices...);
    }

    dtype load(long index) const
    {
      return fast(index);
    }
    // comparison
    template <class K>
    bool operator==(list<K> const &other) const;
    bool operator==(empty_list const &other) const;

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<sliced_list>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif

    // other operations
    template <class V>
    bool contains(V const &v) const;
    intptr_t id() const;

    intptr_t baseid() const
    {
      return reinterpret_cast<intptr_t>(&(*_data));
    }

    long count(T const &x) const;
    template <class Tp, class Sp>
    friend std::ostream &operator<<(std::ostream &os, sliced_list<Tp, Sp> const &v);
  };

  /* list */
  template <class T>
  class list
  {
    static constexpr size_t DEFAULT_CAPACITY = 16;

    // data holder
    using _type = std::remove_cv_t<std::remove_reference_t<T>>;
    typedef container<_type> container_type;
    utils::shared_ref<container_type> _data;

    template <class U, class S>
    friend class sliced_list;

    template <class U>
    friend class list;

  public:
    // types
    typedef typename container_type::value_type value_type;
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::allocator_type allocator_type;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
    typedef typename container_type::reverse_iterator reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;

    // minimal ndarray interface
    typedef typename utils::nested_container_value_type<list>::type dtype;
    static const size_t value = utils::nested_container_depth<list>::value;
    static const bool is_vectorizable = types::is_vectorizable<dtype>::value;
    static const bool is_flat = true;
    static const bool is_strided = false;

    // constructors
    list();
    template <class InputIterator>
    list(InputIterator start, InputIterator stop);
    list(empty_list const &);
    list(size_type sz);
    list(std::initializer_list<T> l);
    list(list<T> &&other);
    list(list<T> const &other);
    template <class F>
    list(list<F> const &other);
    template <class Tp, class S>
    list(sliced_list<Tp, S> const &other);
    template <class Tp, size_t N>
    list(static_list<Tp, N> const &other) : list(other.begin(), other.end())
    {
    }
    template <class Tp, size_t N, class... S>
    list(numpy_gexpr<static_list<Tp, N>, S...> const &other) : list(other.begin(), other.end())
    {
    }
    list<T> &operator=(list<T> &&other);
    template <class F>
    list<T> &operator=(list<F> const &other);
    list<T> &operator=(list<T> const &other);
    list<T> &operator=(empty_list const &);
    template <class Tp, size_t N, class V>
    list<T> &operator=(array_base<Tp, N, V> const &);
    template <class Tp, class S>
    list<T> &operator=(sliced_list<Tp, S> const &other);

    template <class pS>
    list &operator=(ndarray<T, pshape<pS>> const &); // implemented in ndarray.hpp

    template <class S>
    list<T> &operator+=(sliced_list<T, S> const &other);
    template <class S>
    list<T> operator+(sliced_list<T, S> const &other) const;
    template <size_t N, class V>
    list<T> operator+(array_base<T, N, V> const &other) const;

    // io
    template <class S>
    friend std::ostream &operator<<(std::ostream &os, list<S> const &v);

    // comparison
    template <class K>
    bool operator==(list<K> const &other) const;
    bool operator==(empty_list const &) const;
    template <class K>
    bool operator!=(list<K> const &other) const;
    bool operator!=(empty_list const &) const;

    // iterators
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    reverse_iterator rbegin();
    const_reverse_iterator rbegin() const;
    reverse_iterator rend();
    const_reverse_iterator rend() const;

    // comparison
    bool operator<(list<T> const &other) const;
    bool operator<=(list<T> const &other) const;
    bool operator>(list<T> const &other) const;
    bool operator>=(list<T> const &other) const;

// element access
#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<list>;
    using simd_iterator_nobroadcast = simd_iterator;
    template <class vectorizer>
    simd_iterator vbegin(vectorizer) const;
    template <class vectorizer>
    simd_iterator vend(vectorizer) const;
#endif
    reference fast(long n);
    reference operator[](long n);

    const_reference fast(long n) const;
    const_reference operator[](long n) const;

    template <class Sp>
    std::enable_if_t<is_slice<Sp>::value, sliced_list<T, Sp>> operator[](Sp const &s) const;

    template <class... Indices>
    dtype load(long index0, long index1, Indices... indices) const
    {
      return fast(index0).load(index1, indices...);
    }

    dtype load(long index) const
    {
      return fast(index);
    }

    dtype *data()
    {
      return _data->data();
    }
    const dtype *data() const
    {
      return _data->data();
    }

    // modifiers
    template <class Tp>
    void push_back(Tp &&x);
    template <class Tp>
    void insert(long i, Tp &&x);

    void reserve(size_t n);
    void resize(size_t n);
    iterator erase(size_t n);

    T pop(long x = -1);
    void clear();

    // TODO: have to raise a valueError
    none_type remove(T const &x);

    // Misc
    // TODO: have to raise a valueError
    long index(T const &x) const;

    // list interface
    explicit operator bool() const;

    template <class F>
    list<typename __combined<T, F>::type> operator+(list<F> const &s) const;

    template <class F, class S>
    list<decltype(std::declval<T>() + std::declval<typename sliced_list<F, S>::value_type>())>
    operator+(sliced_list<F, S> const &s) const;

    list<T> operator+(empty_list const &) const;
    list<T> operator*(long t) const;
    list<T> const &operator*=(long t);

    template <class F>
    list<T> &operator+=(F const &s);

    long size() const;
    template <class E>
    long _flat_size(E const &e, utils::int_<1>) const;
    template <class E, size_t L>
    long _flat_size(E const &e, utils::int_<L>) const;
    long flat_size() const;

    template <class V>
    bool contains(V const &v) const;
    intptr_t id() const;

    long count(T const &x) const;
    using shape_t = array_tuple<long, value>;
    template <size_t I>
    long shape() const
    {
      if (I == 0)
        return size();
      else
        return details::extract_shape(*this, utils::int_<I>{});
    }

    template <class Tp, size_t N, class V>
    operator array_base<Tp, N, V>() const
    {
      assert(size() == N && "consistent size");
      array_base<Tp, N, V> res;
      std::copy(begin(), end(), res.begin());
      return res;
    }
  };

  template <class T, size_t N>
  list<T> operator*(static_list<T, N> const &self, long t)
  {
    list<T> res(self);
    res *= t;
    return res;
  }
  template <class T, size_t N>
  list<T> operator*(long t, static_list<T, N> const &self)
  {
    return self * t;
  }

  template <class T0, size_t N, class T1>
  list<typename __combined<T0, T1>::type> operator+(static_list<T0, N> const &l0,
                                                    list<T1> const &l1)
  {
    list<typename __combined<T0, T1>::type> out(l0.begin(), l0.end());
    return out += l1;
  }

  /* empty list implementation */
  struct empty_list {
    // minimal ndarray interface
    typedef char dtype;
    static const size_t value = 1;
    static const bool is_vectorizable = false;
    static const bool is_strided = false;
    using shape_t = types::array_tuple<long, value>;
    typedef char value_type;

    typedef empty_iterator iterator;
    typedef empty_iterator const_iterator;
#ifdef USE_XSIMD
    typedef empty_iterator simd_iterator;
    typedef empty_iterator simd_iterator_nobroadcast;
#endif
    template <class T>
    list<T> operator+(list<T> const &s) const;
    template <class T, class S>
    sliced_list<T, S> operator+(sliced_list<T, S> const &s) const;
    template <class T, size_t N, class V>
    static_list<T, N> operator+(array_base<T, N, V> const &s) const;
    empty_list operator+(empty_list const &) const;
    template <class F>
    std::enable_if_t<!is_numexpr_arg<F>::value, list<typename F::value_type>> operator+(F s) const;
    explicit operator bool() const;
    template <class T>
    operator list<T>() const;
    static constexpr long size();

    template <size_t I>
    std::integral_constant<long, 0> shape() const
    {
      return {};
    }

    char fast(long) const
    {
      return {};
    }
    char operator[](long) const
    {
      return {};
    }
    template <class S>
    std::enable_if_t<is_slice<S>::value, empty_list> operator[](S) const
    {
      return {};
    }

    empty_iterator begin() const
    {
      return {};
    }
    empty_iterator end() const
    {
      return {};
    }
  };

  std::ostream &operator<<(std::ostream &os, empty_list const &);
  template <class T, size_t N>
  list<T> operator+(static_list<T, N> const &self, list<T> const &other)
  {
    list<T> res(self.begin(), self.end());
    return res += other;
  }
} // namespace types

namespace utils
{
  /**
   * Reserve enough space to save all values generated from f.
   *
   * We use a dummy arguments (p) to reserve only when f have a
   * const_iterator type.
   */
  template <class T, class From>
  void reserve(types::list<T> &l, From const &f, typename From::const_iterator *p = nullptr);
} // namespace utils

template <class T>
struct assignable<types::list<T>> {
  typedef types::list<typename assignable<T>::type> type;
};

template <class T, class S>
struct assignable<types::sliced_list<T, S>> {
  typedef types::list<typename assignable<T>::type> type;
};

// to cope with std::vector<bool> specialization
template <>
struct returnable<types::list<bool>::reference> {
  using type = bool;
};
PYTHONIC_NS_END

/* overload std::get */
namespace std
{
  template <size_t I, class T>
  typename pythonic::types::list<T>::reference get(pythonic::types::list<T> &t);

  template <size_t I, class T>
  typename pythonic::types::list<T>::const_reference get(pythonic::types::list<T> const &t);

  template <size_t I, class T>
  typename pythonic::types::list<T>::value_type get(pythonic::types::list<T> &&t);

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_list<T, S>::reference get(pythonic::types::sliced_list<T, S> &t);

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_list<T, S>::const_reference
  get(pythonic::types::sliced_list<T, S> const &t);

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_list<T, S>::value_type
  get(pythonic::types::sliced_list<T, S> &&t);

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::list<T>> {
    typedef typename pythonic::types::list<T>::value_type type;
  };
  template <size_t I, class T, class S>
  struct tuple_element<I, pythonic::types::sliced_list<T, S>> {
    typedef typename pythonic::types::sliced_list<T, S>::value_type type;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class A>
struct __combined<container<A>, pythonic::types::empty_list> {
  typedef pythonic::types::list<A> type;
};

template <class A>
struct __combined<pythonic::types::empty_list, container<A>> {
  typedef pythonic::types::list<A> type;
};

template <class A, class B>
struct __combined<container<A>, pythonic::types::list<B>> {
  typedef pythonic::types::list<typename __combined<A, B>::type> type;
};

template <class A, class B>
struct __combined<pythonic::types::list<B>, container<A>> {
  typedef pythonic::types::list<typename __combined<A, B>::type> type;
};

template <class K, class V>
struct __combined<indexable<K>, pythonic::types::list<V>> {
  typedef pythonic::types::list<V> type;
};

template <class V, class K>
struct __combined<pythonic::types::list<V>, indexable<K>> {
  typedef pythonic::types::list<V> type;
};

template <class K, class V0, class V1>
struct __combined<indexable_container<K, V0>, pythonic::types::list<V1>> {
  typedef pythonic::types::list<typename __combined<V0, V1>::type> type;
};

template <class K, class V0, class V1>
struct __combined<pythonic::types::list<V1>, indexable_container<K, V0>> {
  typedef pythonic::types::list<typename __combined<V0, V1>::type> type;
};

template <class K, class V>
struct __combined<indexable_container<K, V>, pythonic::types::empty_list> {
  typedef pythonic::types::list<V> type;
};

template <class K, class V>
struct __combined<pythonic::types::empty_list, indexable_container<K, V>> {
  typedef pythonic::types::list<V> type;
};

template <class T0, class T1>
struct __combined<pythonic::types::list<T0>, pythonic::types::list<T1>> {
  typedef pythonic::types::list<typename __combined<T0, T1>::type> type;
};

template <class T, class S>
struct __combined<pythonic::types::sliced_list<T, S>, pythonic::types::empty_list> {
  typedef pythonic::types::list<T> type;
};
template <class T, class S>
struct __combined<pythonic::types::empty_list, pythonic::types::sliced_list<T, S>> {
  typedef pythonic::types::list<T> type;
};

template <class T0, class T1, class S>
struct __combined<pythonic::types::sliced_list<T1, S>, pythonic::types::list<T0>> {
  typedef pythonic::types::list<typename __combined<T0, T1>::type> type;
};
template <class T0, class T1, class S>
struct __combined<pythonic::types::list<T0>, pythonic::types::sliced_list<T1, S>> {
  typedef pythonic::types::list<typename __combined<T0, T1>::type> type;
};

template <class T, size_t N, class V>
struct __combined<pythonic::types::array_base<T, N, V>, pythonic::types::empty_list> {
  typedef pythonic::types::list<T> type;
};
template <class T, size_t N, class V>
struct __combined<pythonic::types::empty_list, pythonic::types::array_base<T, N, V>> {
  typedef pythonic::types::list<T> type;
};
template <class T, size_t N, class V, class Tp>
struct __combined<pythonic::types::array_base<T, N, V>, pythonic::types::list<Tp>> {
  typedef pythonic::types::list<typename __combined<T, Tp>::type> type;
};
template <class T, size_t N, class V, class Tp>
struct __combined<pythonic::types::list<Tp>, pythonic::types::array_base<T, N, V>> {
  typedef pythonic::types::list<typename __combined<T, Tp>::type> type;
};

/* } */

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

template <>
struct to_python<typename std::vector<bool>::reference> {
  static PyObject *convert(typename std::vector<bool>::reference const &v);
};

struct phantom_type; // ghost don't exist
template <>
struct to_python<
    std::conditional_t<std::is_same<bool, typename std::vector<bool>::const_reference>::value,
                       phantom_type, typename std::vector<bool>::const_reference>> {
  static PyObject *convert(typename std::vector<bool>::const_reference const &v);
};

template <typename T>
struct to_python<types::list<T>> {
  static PyObject *convert(types::list<T> const &v);
};
template <typename T, typename S>
struct to_python<types::sliced_list<T, S>> {
  static PyObject *convert(types::sliced_list<T, S> const &v);
};
template <>
struct to_python<types::empty_list> {
  static PyObject *convert(types::empty_list const &);
};

template <class T>
struct from_python<types::list<T>> {

  static bool is_convertible(PyObject *obj);

  static types::list<T> convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif

#endif
