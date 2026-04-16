#ifndef PYTHONIC_TYPES_ARRAY_HPP
#define PYTHONIC_TYPES_ARRAY_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/types/nditerator.hpp"

#include "pythonic/builtins/len.hpp"
#include "pythonic/types/bool.hpp"
#include "pythonic/types/slice.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/reserve.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include "pythonic/builtins/NotImplementedError.hpp"

#include <algorithm>
#include <cassert>

PYTHONIC_NS_BEGIN

namespace types
{

  /// Sliced array

  // Constructors
  template <class T, class S>
  sliced_array<T, S>::sliced_array() : _data(utils::no_memory())
  {
  }
  template <class T, class S>
  sliced_array<T, S>::sliced_array(sliced_array<T, S> const &s) : _data(s._data), slicing(s.slicing)
  {
  }
  template <class T, class S>
  sliced_array<T, S>::sliced_array(array<T> const &other, S const &s)
      : _data(other._data), slicing(s.normalize(other.size()))
  {
  }
  template <class T, class S>
  template <class Sn>
  sliced_array<T, S>::sliced_array(utils::shared_ref<container_type> const &other, Sn const &s)
      : _data(other), slicing(s)
  {
  }

  // iterators
  template <class T, class S>
  typename sliced_array<T, S>::iterator sliced_array<T, S>::begin()
  {
    return {*this, 0};
  }
  template <class T, class S>
  typename sliced_array<T, S>::const_iterator sliced_array<T, S>::begin() const
  {
    return {*this, 0};
  }
  template <class T, class S>
  typename sliced_array<T, S>::iterator sliced_array<T, S>::end()
  {
    return {*this, size()};
  }
  template <class T, class S>
  typename sliced_array<T, S>::const_iterator sliced_array<T, S>::end() const
  {
    return {*this, size()};
  }

  // size
  template <class T, class S>
  long sliced_array<T, S>::size() const
  {
    return slicing.size();
  }
  template <class T, class S>
  sliced_array<T, S>::operator bool() const
  {
    return slicing.size();
  }

  // accessor
  template <class T, class S>
  typename sliced_array<T, S>::const_reference sliced_array<T, S>::fast(long i) const
  {
    assert(0 <= i && i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }
  template <class T, class S>
  typename sliced_array<T, S>::reference sliced_array<T, S>::fast(long i)
  {
    assert(0 <= i && i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }
  template <class T, class S>
  typename sliced_array<T, S>::const_reference sliced_array<T, S>::operator[](long i) const
  {
    assert(i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }
  template <class T, class S>
  typename sliced_array<T, S>::reference sliced_array<T, S>::operator[](long i)
  {
    assert(i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }

  template <class T, class S>
  template <class Sp>
  std::enable_if_t<is_slice<Sp>::value,
                   sliced_array<T, decltype(std::declval<S>() * std::declval<Sp>())>>
  sliced_array<T, S>::operator[](Sp s) const
  {
    return {_data, slicing * s.normalize(this->size())};
  }

  // io
  template <class Tp, class Sp>
  std::ostream &operator<<(std::ostream &os, sliced_array<Tp, Sp> const &v)
  {
    os << '[';
    auto iter = v.begin();
    if (iter != v.end()) {
      while (iter + 1 != v.end()) {
        os << *iter << ", ";
        ++iter;
      }
      os << *iter;
    }
    return os << ']';
  }

  // comparison
  template <class T, class S>
  template <class K>
  bool sliced_array<T, S>::operator==(array<K> const &other) const
  {
    if (size() != other.size())
      return false;
    return std::equal(begin(), end(), other.begin());
  }
  template <class T, class S>
  inline sliced_array<T, S> &sliced_array<T, S>::operator=(sliced_array<T, S> const &s)
  {
    if (slicing.step == 1) {
      // inserting before erasing in case of self-copy
      auto insert_pt = _data->begin() + slicing.lower;
      _data->insert(insert_pt, s.begin(), s.end());
      auto erase_pt = _data->begin() + s.size();
      _data->erase(erase_pt + slicing.lower, erase_pt + slicing.upper);
    } else
      assert(!"not implemented yet");
    return *this;
  }
  template <class T, class S>
  sliced_array<T, S> &sliced_array<T, S>::operator=(array<T> const &seq)
  {
    if (slicing.step == 1) {
      // inserting before erasing in case of self-copy
      auto insert_pt = _data->begin() + slicing.lower;
      _data->insert(insert_pt, seq.begin(), seq.end());
      auto erase_pt = _data->begin() + seq.size();
      _data->erase(erase_pt + slicing.lower, erase_pt + slicing.upper);
    } else
      assert(!"not implemented yet");
    return *this;
  }
  template <class T, class S>
  array<T> sliced_array<T, S>::operator+(array<T> const &s) const
  {
    array<T> out(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), out.begin()));
    return out;
  }
  template <class T, class S>
  template <size_t N, class V>
  array<T> sliced_array<T, S>::operator+(array_base<T, N, V> const &s) const
  {
    array<T> out(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), out.begin()));
    return out;
  }
  template <class T, class S>
  template <class Tp, class Sp>
  array<typename __combined<T, Tp>::type>
  sliced_array<T, S>::operator+(sliced_array<Tp, Sp> const &s) const
  {
    array<typename __combined<T, Tp>::type> out(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), out.begin()));
    return out;
  }
  template <class N, class T>
  array<T> operator*(N n, array<T> const &l)
  {
    return l * n;
  }
#ifdef USE_XSIMD
  template <class T, class S>
  template <class vectorizer>
  typename sliced_array<T, S>::simd_iterator sliced_array<T, S>::vbegin(vectorizer) const
  {
    return {_data->data() + slicing.lower};
  }

  template <class T, class S>
  template <class vectorizer>
  typename sliced_array<T, S>::simd_iterator sliced_array<T, S>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {_data->data() + slicing.lower + long(size() / vector_size * vector_size)};
  }

#endif

  // other operations
  template <class T, class S>
  template <class V>
  bool sliced_array<T, S>::contains(V const &v) const
  {
    return std::find(_data->begin(), _data->end(), v) != _data->end();
  }
  template <class T, class S>
  intptr_t sliced_array<T, S>::id() const
  {
    // sharing is not implemented for sliced array
    return reinterpret_cast<intptr_t>(this);
  }

  template <class T, class S>
  long sliced_array<T, S>::count(T const &x) const
  {
    return std::count(begin(), end(), x);
  }

  /// List

  // constructors
  template <class T>
  array<T>::array() : _data(utils::no_memory())
  {
  }
  template <class T>
  template <class InputIterator>
  array<T>::array(InputIterator start, InputIterator stop) : _data()
  {
    if (std::is_same<typename std::iterator_traits<InputIterator>::iterator_category,
                     std::random_access_iterator_tag>::value)
      _data->reserve(std::distance(start, stop));
    else
      _data->reserve(DEFAULT_CAPACITY);
    std::copy(start, stop, std::back_inserter(*_data));
  }
  template <class T>
  array<T>::array(size_type sz) : _data(sz)
  {
  }
  template <class T>
  array<T>::array(array<T> &&other) : _data(std::move(other._data))
  {
  }
  template <class T>
  array<T>::array(array<T> const &other) : _data(other._data)
  {
  }
  template <class T>
  template <class F>
  array<T>::array(array<F> const &other) : _data(other.size())
  {
    std::copy(other.begin(), other.end(), begin());
  }
  template <class T>
  template <class Tp, class S>
  array<T>::array(sliced_array<Tp, S> const &other) : _data(other.begin(), other.end())
  {
  }

  // operators
  template <class T>
  array<T> &array<T>::operator=(array<T> &&other)
  {
    _data = std::move(other._data);
    return *this;
  }
  template <class T>
  template <class F>
  array<T> &array<T>::operator=(array<F> const &other)
  {
    _data = utils::shared_ref<container_type>{other.size()};
    std::copy(other.begin(), other.end(), begin());
    return *this;
  }
  template <class T>
  array<T> &array<T>::operator=(array<T> const &other)
  {
    _data = other._data;
    return *this;
  }
  template <class T>
  template <class Tp, size_t N, class V>
  array<T> &array<T>::operator=(array_base<Tp, N, V> const &other)
  {
    _data = utils::shared_ref<container_type>(other.begin(), other.end());
    return *this;
  }
  template <class T>
  template <class Tp, class S>
  array<T> &array<T>::operator=(sliced_array<Tp, S> const &other)
  {
    if (other._data == _data) {
      auto it = std::copy(other.begin(), other.end(), _data->begin());
      _data->resize(it - _data->begin());
    } else
      _data = utils::shared_ref<container_type>(other.begin(), other.end());
    return *this;
  }

  template <class T>
  template <class S>
  array<T> &array<T>::operator+=(sliced_array<T, S> const &other)
  {
    _data->resize(size() + other.size());
    std::copy(other.begin(), other.end(), _data->begin());
    return *this;
  }

  template <class T>
  template <class S>
  array<T> array<T>::operator+(sliced_array<T, S> const &other) const
  {
    array<T> new_array(begin(), end());
    new_array.reserve(size() + other.size());
    std::copy(other.begin(), other.end(), std::back_inserter(new_array));
    return new_array;
  }

  template <class T>
  template <size_t N, class V>
  array<T> array<T>::operator+(array_base<T, N, V> const &other) const
  {
    array<T> new_array(begin(), end());
    new_array.reserve(size() + other.size());
    std::copy(other.begin(), other.end(), std::back_inserter(new_array));
    return new_array;
  }

  // io
  template <class T>
  std::ostream &operator<<(std::ostream &os, array<T> const &v)
  {
    os << '[';
    auto iter = v.begin();
    if (iter != v.end()) {
      while (iter + 1 != v.end())
        os << *iter++ << ", ";
      os << *iter;
    }
    return os << ']';
  }

  // comparison
  template <class T>
  template <class K>
  bool array<T>::operator==(array<K> const &other) const
  {
    if (size() != other.size())
      return false;
    return std::equal(begin(), end(), other.begin());
  }
  template <class T>
  template <class K>
  bool array<T>::operator!=(array<K> const &other) const
  {
    return !operator==(other);
  }

  // iterators
  template <class T>
  typename array<T>::iterator array<T>::begin()
  {
    return {*this, 0};
  }
  template <class T>
  typename array<T>::const_iterator array<T>::begin() const
  {
    return {*this, 0};
  }
  template <class T>
  typename array<T>::iterator array<T>::end()
  {
    return {*this, size()};
  }
  template <class T>
  typename array<T>::const_iterator array<T>::end() const
  {
    return {*this, size()};
  }

  // comparison
  template <class T>
  bool array<T>::operator<(array<T> const &other) const
  {
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end());
  }
  template <class T>
  bool array<T>::operator>(array<T> const &other) const
  {
    return std::lexicographical_compare(other.begin(), other.end(), begin(), end());
  }
  template <class T>
  bool array<T>::operator<=(array<T> const &other) const
  {
    return !(*this > other);
  }
  template <class T>
  bool array<T>::operator>=(array<T> const &other) const
  {
    return !(*this < other);
  }

// element access
#ifdef USE_XSIMD
  template <class T>
  template <class vectorizer>
  typename array<T>::simd_iterator array<T>::vbegin(vectorizer) const
  {
    return {_data->data()};
  }

  template <class T>
  template <class vectorizer>
  typename array<T>::simd_iterator array<T>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {_data->data() + long(size() / vector_size * vector_size)};
  }

#endif
  template <class T>
  typename array<T>::reference array<T>::fast(long n)
  {
    return (*_data)[n];
  }
  template <class T>
  typename array<T>::reference array<T>::operator[](long n)
  {
    if (n < 0)
      n += size();
    assert(0 <= n && n < size());
    return fast(n);
  }
  template <class T>
  typename array<T>::const_reference array<T>::fast(long n) const
  {
    assert(n < size());
    return (*_data)[n];
  }
  template <class T>
  typename array<T>::const_reference array<T>::operator[](long n) const
  {
    if (n < 0)
      n += size();
    assert(0 <= n && n < size());
    return fast(n);
  }

  template <class T>
  template <class Sp>
  std::enable_if_t<is_slice<Sp>::value, sliced_array<T, Sp>> array<T>::operator[](Sp const &s) const
  {
    return {*this, s};
  }

  // modifiers
  template <class T>
  template <class Tp>
  void array<T>::push_back(Tp &&x)
  {
    // FIXME: clang-3.4 doesn't support emplace_back for vector of bool
    _data->push_back(std::forward<Tp>(x));
  }
  template <class T>
  template <class Tp>
  void array<T>::insert(long i, Tp &&x)
  {
    if (i == size())
      _data->emplace_back(std::forward<Tp>(x));
    else
      _data->insert(_data->begin() + i, std::forward<Tp>(x));
  }
  template <class T>
  void array<T>::reserve(size_t n)
  {
    if (n > _data->capacity())
      _data->reserve((n / 2) * 3);
  }
  template <class T>
  void array<T>::resize(size_t n)
  {
    _data->resize(n);
  }
  template <class T>
  void array<T>::erase(size_t n)
  {
    _data->erase(_data->begin() + n);
  }
  template <class T>
  T array<T>::pop(long x)
  {
    long sz = size();
    x = x % sz;
    if (x < 0)
      x += sz;
    T res = fast(x);
    erase(x);
    return res;
  }
  template <class T>
  void array<T>::clear()
  {
    _data->clear();
  }

  // TODO: have to raise a valueError
  template <class T>
  none_type array<T>::remove(T const &x)
  {
    erase(index(x));
    return {};
  }

  // Misc
  template <class T>
  long array<T>::index(T const &x) const
  {
    return std::find(begin(), end(), x) - begin();
  }

  // array interface
  template <class T>
  array<T>::operator bool() const
  {
    return !_data->empty();
  }

  template <class T>
  template <class F>
  array<typename __combined<T, F>::type> array<T>::operator+(array<F> const &s) const
  {
    array<typename __combined<T, F>::type> clone(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), clone.begin()));
    return clone;
  }

  template <class T>
  template <class F, class S>
  array<decltype(std::declval<T>() + std::declval<typename sliced_array<F, S>::value_type>())>
  array<T>::operator+(sliced_array<F, S> const &s) const
  {
    array<decltype(std::declval<T>() + std::declval<typename sliced_array<F, S>::value_type>())>
        clone(size() + len(s));
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), clone.begin()));
    return clone;
  }

  template <class T>
  array<T> array<T>::operator*(long n) const
  {
    if (size() == 1) {
      array<T> r(n);
      std::fill(begin(), end(), fast(0));
      return r;
    } else {
      array<T> r(size() * n);
      auto start = r.begin();
      while (start != r.end())
        start = std::copy(this->begin(), this->end(), start);
      return r;
    }
  }
  template <class T>
  array<T> const &array<T>::operator*=(long n)
  {
    if (size() == 1) {
      resize(n);
      std::fill(begin() + 1, end(), fast(0));
    } else {
      auto const initial_size = size();
      resize(n * initial_size);
      // FIXME: could use less calls to std::copy
      auto tgt = begin() + initial_size;
      for (long i = 1; i < n; ++i)
        tgt = std::copy(begin(), begin() + initial_size, tgt);
    }
    return *this;
  }

  template <class T>
  template <class F>
  array<T> &array<T>::operator+=(F const &s)
  {
    reserve(size() + s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(*this));
    return *this;
  }

  template <class T>
  long array<T>::size() const
  {
    return _data->size();
  }
  template <class T>
  template <class E>
  long array<T>::_flat_size(E const &e, utils::int_<1>) const
  {
    return std::distance(e.begin(), e.end());
  }
  template <class T>
  template <class E, size_t L>
  long array<T>::_flat_size(E const &e, utils::int_<L>) const
  {
    return std::distance(e.begin(), e.end()) * _flat_size(e[0], utils::int_<L - 1>{});
  }
  template <class T>
  long array<T>::flat_size() const
  {
    return _flat_size(*this, utils::int_<value>{});
  }

  template <class T>
  template <class V>
  bool array<T>::contains(V const &v) const
  {
    return std::find(_data->begin(), _data->end(), v) != _data->end();
  }
  template <class T>
  intptr_t array<T>::id() const
  {
    return reinterpret_cast<intptr_t>(&(*_data));
  }

  template <class T>
  long array<T>::count(T const &x) const
  {
    return std::count(begin(), end(), x);
  }

} // namespace types

namespace utils
{

  template <class T, class From>
  void reserve(types::array<T> &l, From const &f, typename From::const_iterator *)
  {
    l.reserve(builtins::len(f));
  }
} // namespace utils
PYTHONIC_NS_END

/* overload std::get */
namespace std
{
  template <size_t I, class T>
  typename pythonic::types::array<T>::reference get(pythonic::types::array<T> &t)
  {
    return t[I];
  }

  template <size_t I, class T>
  typename pythonic::types::array<T>::const_reference get(pythonic::types::array<T> const &t)
  {
    return t[I];
  }

  template <size_t I, class T>
  typename pythonic::types::array<T>::value_type get(pythonic::types::array<T> &&t)
  {
    return std::move(t)[I];
  }

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_array<T, S>::reference
  get(pythonic::types::sliced_array<T, S> &t)
  {
    return t[I];
  }

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_array<T, S>::const_reference
  get(pythonic::types::sliced_array<T, S> const &t)
  {
    return t[I];
  }

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_array<T, S>::value_type
  get(pythonic::types::sliced_array<T, S> &&t)
  {
    return std::move(t)[I];
  }
} // namespace std

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

template <class T>
PyObject *to_python<types::array<T>>::convert(types::array<T> const &v)
{
  throw types::NotImplementedError("Pythran cannot efficiently convert array::array values");
}
template <class T, class S>
PyObject *to_python<types::sliced_array<T, S>>::convert(types::sliced_array<T, S> const &v)
{
  throw types::NotImplementedError("Pythran cannot efficiently convert array::array values");
}

PYTHONIC_NS_END

#endif

#endif
