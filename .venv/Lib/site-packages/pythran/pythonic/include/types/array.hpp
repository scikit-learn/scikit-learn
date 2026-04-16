#ifndef PYTHONIC_INCLUDE_TYPES_ARRAY_HPP
#define PYTHONIC_INCLUDE_TYPES_ARRAY_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/empty_iterator.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/nditerator.hpp"
#include "pythonic/include/types/slice.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/types/vectorizable_type.hpp"
#include "pythonic/include/utils/allocate.hpp"
#include "pythonic/include/utils/int_.hpp"
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
  template <class T>
  class array;
  template <class T, class S>
  class sliced_array;
  template <class T, class pS>
  struct ndarray;
  template <class... Tys>
  struct pshape;

  template <class E>
  struct array_reference {
    typename E::data_type &data;
    array_reference(typename E::data_type &data) : data(data)
    {
    }
    array_reference &operator=(typename E::value_type value)
    {
      data = value;
      return *this;
    }
    array_reference &operator=(array_reference value)
    {
      data = value.data;
      return *this;
    }

    operator typename E::value_type() const
    {
      return data;
    }

    friend void swap(array_reference self, array_reference other)
    {
      std::swap(self.data, other.data);
    }
  };

  template <class E>
  struct array_iterator : std::iterator<std::random_access_iterator_tag, typename E::value_type> {
    E data;
    long index;
    array_iterator(E data, long index) : data(data), index(index)
    {
    }

    array_reference<E> operator*()
    {
      return data.fast(index);
    }
    auto operator*() const -> decltype(data.fast(index))
    {
      return data.fast(index);
    }
    array_iterator &operator++()
    {
      ++index;
      return *this;
    }
    array_iterator &operator--()
    {
      --index;
      return *this;
    }
    array_iterator &operator+=(long i)
    {
      index += i;
      return *this;
    }
    array_iterator &operator-=(long i)
    {
      index -= i;
      return *this;
    }
    array_iterator operator+(long i) const
    {
      array_iterator res(*this);
      res += i;
      return res;
    }
    array_iterator operator-(long i) const
    {
      array_iterator res(*this);
      res -= i;
      return res;
    }
    long operator-(array_iterator const &other) const
    {
      return index - other.index;
    }
    bool operator!=(array_iterator const &other) const
    {
      return index != other.index;
    }
    bool operator==(array_iterator const &other) const
    {
      return index == other.index;
    }
    bool operator<(array_iterator const &other) const
    {
      return index < other.index;
    }
    bool operator<=(array_iterator const &other) const
    {
      return index <= other.index;
    }
    array_iterator &operator=(array_iterator const &other)
    {
      index = other.index;
      return *this;
    }
  };

  /* array view */
  template <class T, class S = slice>
  class sliced_array
  {

    // data holder
    using _type = std::remove_cv_t<std::remove_reference_t<T>>;
    typedef container<_type> container_type;
    utils::shared_ref<container_type> _data;

    template <class U>
    friend class array;

    typename S::normalized_type slicing;

  public:
    //  types
    typedef T data_type;
    typedef std::conditional_t<std::is_integral<T>::value, long, double> value_type;
    typedef array_reference<sliced_array> reference;
    typedef value_type const_reference;
    typedef array_iterator<sliced_array> iterator;
    typedef array_iterator<sliced_array> const_iterator;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::allocator_type allocator_type;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    // minimal ndarray interface
    typedef data_type dtype;
    static const size_t value = 1;
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
    sliced_array();
    sliced_array(sliced_array<T, S> const &s);
    sliced_array(array<T> const &other, S const &s);
    template <class Sn>
    sliced_array(utils::shared_ref<container_type> const &other, Sn const &s);

    // assignment
    sliced_array &operator=(array<T> const &);
    sliced_array &operator=(sliced_array<T, S> const &);
    array<T> operator+(array<T> const &) const;
    template <size_t N, class V>
    array<T> operator+(array_base<T, N, V> const &) const;
    template <class Tp, class Sp>
    array<typename __combined<T, Tp>::type> operator+(sliced_array<Tp, Sp> const &) const;

    // iterators
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    reverse_iterator rbegin()
    {
      return {end()};
    }
    const_reverse_iterator rbegin() const
    {
      return {end()};
    }
    reverse_iterator rend()
    {
      return {begin()};
    }
    const_reverse_iterator rend() const
    {
      return {begin()};
    }

    // size
    long size() const;
    explicit operator bool() const;

    // accessors
    const_reference fast(long i) const;
    reference fast(long i);
    const_reference operator[](long i) const;
    reference operator[](long i);
    template <class Sp>
    std::enable_if_t<is_slice<Sp>::value,
                     sliced_array<T, decltype(std::declval<S>() * std::declval<Sp>())>>
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
    bool operator==(array<K> const &other) const;

#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<sliced_array>;
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
    friend std::ostream &operator<<(std::ostream &os, sliced_array<Tp, Sp> const &v);
  };

  /* array */
  template <class T>
  class array
  {

    static const size_t DEFAULT_CAPACITY = 16;

    // data holder
    using _type = std::remove_cv_t<std::remove_reference_t<T>>;
    typedef container<_type> container_type;
    utils::shared_ref<container_type> _data;

    template <class U, class S>
    friend class sliced_array;

    template <class U>
    friend class array;

  public:
    // types
    typedef T data_type;
    typedef std::conditional_t<std::is_integral<T>::value, long, double> value_type;
    typedef array_reference<array> reference;
    typedef value_type const_reference;
    typedef array_iterator<array> iterator;
    typedef array_iterator<array> const_iterator;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::allocator_type allocator_type;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    // minimal ndarray interface
    typedef data_type dtype;
    static const size_t value = 1;
    static const bool is_vectorizable = types::is_vectorizable<dtype>::value;
    static const bool is_flat = true;
    static const bool is_strided = false;

    // constructors
    array();
    template <class InputIterator>
    array(InputIterator start, InputIterator stop);
    array(size_type sz);
    array(array &&other);
    array(array const &other);
    template <class Tp>
    array(list<Tp> const &other) : array(other.begin(), other.end())
    {
    }
    template <class Tp, size_t N>
    array(static_list<Tp, N> const &other) : array(other.begin(), other.end())
    {
    }
    template <class F>
    array(array<F> const &other);
    template <class Tp, class S>
    array(sliced_array<Tp, S> const &other);
    array<T> &operator=(array<T> &&other);
    template <class F>
    array<T> &operator=(array<F> const &other);
    array<T> &operator=(array<T> const &other);
    template <class Tp, size_t N, class V>
    array<T> &operator=(array_base<Tp, N, V> const &);
    template <class Tp, class S>
    array<T> &operator=(sliced_array<Tp, S> const &other);

    template <class pS>
    array &operator=(ndarray<T, pshape<pS>> const &); // implemented in ndarray.hpp

    template <class S>
    array<T> &operator+=(sliced_array<T, S> const &other);
    template <class S>
    array<T> operator+(sliced_array<T, S> const &other) const;
    template <size_t N, class V>
    array<T> operator+(array_base<T, N, V> const &other) const;

    // io
    template <class S>
    friend std::ostream &operator<<(std::ostream &os, array<S> const &v);

    // comparison
    template <class K>
    bool operator==(array<K> const &other) const;
    template <class K>
    bool operator!=(array<K> const &other) const;

    // iterators
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    reverse_iterator rbegin()
    {
      return {end()};
    }
    const_reverse_iterator rbegin() const
    {
      return {end()};
    }
    reverse_iterator rend()
    {
      return {begin()};
    }
    const_reverse_iterator rend() const
    {
      return {begin()};
    }

    // comparison
    bool operator<(array<T> const &other) const;
    bool operator<=(array<T> const &other) const;
    bool operator>(array<T> const &other) const;
    bool operator>=(array<T> const &other) const;

// element access
#ifdef USE_XSIMD
    using simd_iterator = const_simd_nditerator<array>;
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
    std::enable_if_t<is_slice<Sp>::value, sliced_array<T, Sp>> operator[](Sp const &s) const;

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
    void erase(size_t n);

    T pop(long x = -1);
    void clear();

    // TODO: have to raise a valueError
    none_type remove(T const &x);

    // Misc
    // TODO: have to raise a valueError
    long index(T const &x) const;

    // array interface
    explicit operator bool() const;

    template <class F>
    array<typename __combined<T, F>::type> operator+(array<F> const &s) const;

    template <class F, class S>
    array<decltype(std::declval<T>() + std::declval<typename sliced_array<F, S>::value_type>())>
    operator+(sliced_array<F, S> const &s) const;

    array<T> operator*(long t) const;
    array<T> const &operator*=(long t);

    template <class F>
    array<T> &operator+=(F const &s);

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
  void reserve(types::array<T> &l, From const &f, typename From::const_iterator *p = nullptr);
} // namespace utils

template <class T>
struct assignable<types::array<T>> {
  typedef types::array<typename assignable<T>::type> type;
};

template <class T, class S>
struct assignable<types::sliced_array<T, S>> {
  typedef types::array<typename assignable<T>::type> type;
};

template <class E>
struct assignable<types::array_reference<E>> {
  typedef typename E::value_type type;
};

PYTHONIC_NS_END

/* overload std::get */
namespace std
{
  template <size_t I, class T>
  typename pythonic::types::array<T>::reference get(pythonic::types::array<T> &t);

  template <size_t I, class T>
  typename pythonic::types::array<T>::const_reference get(pythonic::types::array<T> const &t);

  template <size_t I, class T>
  typename pythonic::types::array<T>::value_type get(pythonic::types::array<T> &&t);

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_array<T, S>::reference
  get(pythonic::types::sliced_array<T, S> &t);

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_array<T, S>::const_reference
  get(pythonic::types::sliced_array<T, S> const &t);

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_array<T, S>::value_type
  get(pythonic::types::sliced_array<T, S> &&t);

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::array<T>> {
    typedef typename pythonic::types::array<T>::value_type type;
  };
  template <size_t I, class T, class S>
  struct tuple_element<I, pythonic::types::sliced_array<T, S>> {
    typedef typename pythonic::types::sliced_array<T, S>::value_type type;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class A, class B>
struct __combined<container<A>, pythonic::types::array<B>> {
  typedef pythonic::types::array<typename __combined<A, B>::type> type;
};

template <class A, class B>
struct __combined<pythonic::types::array<B>, container<A>> {
  typedef pythonic::types::array<typename __combined<A, B>::type> type;
};

template <class K, class V>
struct __combined<indexable<K>, pythonic::types::array<V>> {
  typedef pythonic::types::array<V> type;
};

template <class V, class K>
struct __combined<pythonic::types::array<V>, indexable<K>> {
  typedef pythonic::types::array<V> type;
};

template <class K, class V0, class V1>
struct __combined<indexable_container<K, V0>, pythonic::types::array<V1>> {
  typedef pythonic::types::array<typename __combined<V0, V1>::type> type;
};

template <class K, class V0, class V1>
struct __combined<pythonic::types::array<V1>, indexable_container<K, V0>> {
  typedef pythonic::types::array<typename __combined<V0, V1>::type> type;
};

template <class T0, class T1>
struct __combined<pythonic::types::array<T0>, pythonic::types::array<T1>> {
  typedef pythonic::types::array<typename __combined<T0, T1>::type> type;
};

template <class T0, class T1, class S>
struct __combined<pythonic::types::sliced_array<T1, S>, pythonic::types::array<T0>> {
  typedef pythonic::types::array<typename __combined<T0, T1>::type> type;
};
template <class T0, class T1, class S>
struct __combined<pythonic::types::array<T0>, pythonic::types::sliced_array<T1, S>> {
  typedef pythonic::types::array<typename __combined<T0, T1>::type> type;
};

template <class T, size_t N, class V, class Tp>
struct __combined<pythonic::types::array_base<T, N, V>, pythonic::types::array<Tp>> {
  typedef pythonic::types::array<typename __combined<T, Tp>::type> type;
};
template <class T, size_t N, class V, class Tp>
struct __combined<pythonic::types::array<Tp>, pythonic::types::array_base<T, N, V>> {
  typedef pythonic::types::array<typename __combined<T, Tp>::type> type;
};

/* } */

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

template <typename T>
struct to_python<types::array<T>> {
  static PyObject *convert(types::array<T> const &v);
};
template <typename T, typename S>
struct to_python<types::sliced_array<T, S>> {
  static PyObject *convert(types::sliced_array<T, S> const &v);
};

PYTHONIC_NS_END

#endif

#endif
