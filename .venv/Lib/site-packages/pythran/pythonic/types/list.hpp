#ifndef PYTHONIC_TYPES_LIST_HPP
#define PYTHONIC_TYPES_LIST_HPP

#include "pythonic/include/types/list.hpp"
#include "pythonic/types/nditerator.hpp"

#include "pythonic/builtins/len.hpp"
#include "pythonic/types/bool.hpp"
#include "pythonic/types/slice.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/reserve.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include <algorithm>
#include <cassert>

PYTHONIC_NS_BEGIN

namespace types
{

  /// Sliced list

  // Constructors
  template <class T, class S>
  sliced_list<T, S>::sliced_list() : _data(utils::no_memory())
  {
  }
  template <class T, class S>
  sliced_list<T, S>::sliced_list(sliced_list<T, S> const &s) : _data(s._data), slicing(s.slicing)
  {
  }
  template <class T, class S>
  sliced_list<T, S>::sliced_list(list<T> const &other, S const &s)
      : _data(other._data), slicing(s.normalize(other.size()))
  {
  }
  template <class T, class S>
  template <class Sn>
  sliced_list<T, S>::sliced_list(utils::shared_ref<container_type> const &other, Sn const &s)
      : _data(other), slicing(s)
  {
  }

  // iterators
  template <class T, class S>
  typename sliced_list<T, S>::iterator sliced_list<T, S>::begin()
  {
    return {*this, 0};
  }
  template <class T, class S>
  typename sliced_list<T, S>::const_iterator sliced_list<T, S>::begin() const
  {
    return {*this, 0};
  }
  template <class T, class S>
  typename sliced_list<T, S>::iterator sliced_list<T, S>::end()
  {
    return {*this, size()};
  }
  template <class T, class S>
  typename sliced_list<T, S>::const_iterator sliced_list<T, S>::end() const
  {
    return {*this, size()};
  }

  // size
  template <class T, class S>
  long sliced_list<T, S>::size() const
  {
    return slicing.size();
  }
  template <class T, class S>
  sliced_list<T, S>::operator bool() const
  {
    return slicing.size();
  }

  // accessor
  template <class T, class S>
  typename sliced_list<T, S>::const_reference sliced_list<T, S>::fast(long i) const
  {
    assert(0 <= i && i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }
  template <class T, class S>
  typename sliced_list<T, S>::const_reference sliced_list<T, S>::operator[](long i) const
  {
    assert(i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }
  template <class T, class S>
  typename sliced_list<T, S>::reference sliced_list<T, S>::operator[](long i)
  {
    assert(i < size());
    auto const index = slicing.get(i);
    assert(0 <= index && index < (long)_data->size());
    return (*_data)[index];
  }

  template <class T, class S>
  template <class Sp>
  std::enable_if_t<is_slice<Sp>::value,
                   sliced_list<T, decltype(std::declval<S>() * std::declval<Sp>())>>
  sliced_list<T, S>::operator[](Sp s) const
  {
    return {_data, slicing * s.normalize(this->size())};
  }

  // io
  template <class Tp, class Sp>
  std::ostream &operator<<(std::ostream &os, sliced_list<Tp, Sp> const &v)
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
  bool sliced_list<T, S>::operator==(list<K> const &other) const
  {
    if (size() != other.size())
      return false;
    return std::equal(begin(), end(), other.begin());
  }
  template <class T, class S>
  bool sliced_list<T, S>::operator==(empty_list const &other) const
  {
    return size() == 0;
  }
  template <class T, class S>
  inline sliced_list<T, S> &sliced_list<T, S>::operator=(sliced_list<T, S> const &s)
  {
    if (slicing.step == 1) {
      // no sharing
      if (_data->data() != s._data->data()) {
        auto insert_pt = _data->begin() + slicing.lower;
        _data->insert(insert_pt, s.begin(), s.end());
      }
      // sharing
      else {
        std::vector<T> tmp{s.begin(), s.end()};
        auto insert_pt = _data->begin() + slicing.lower;
        _data->insert(insert_pt, tmp.begin(), tmp.end());
      }
      auto erase_pt = _data->begin() + s.size();
      _data->erase(erase_pt + slicing.lower, erase_pt + slicing.upper);
    } else
      assert(!"not implemented yet");
    return *this;
  }
  template <class T, class S>
  sliced_list<T, S> &sliced_list<T, S>::operator=(list<T> const &seq)
  {
    if (slicing.step == 1) {
      // no sharing
      if (_data->data() != seq._data->data()) {
        auto insert_pt = _data->begin() + slicing.lower;
        _data->insert(insert_pt, seq.begin(), seq.end());
      }
      // sharing
      else {
        std::vector<T> tmp{seq.begin(), seq.end()};
        auto insert_pt = _data->begin() + slicing.lower;
        _data->insert(insert_pt, tmp.begin(), tmp.end());
      }
      auto erase_pt = _data->begin() + seq.size();
      _data->erase(erase_pt + slicing.lower, erase_pt + slicing.upper);
    } else
      assert(!"not implemented yet");
    return *this;
  }
  template <class T, class S>
  list<T> sliced_list<T, S>::operator+(list<T> const &s) const
  {
    list<T> out(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), out.begin()));
    return out;
  }
  template <class T, class S>
  template <size_t N, class V>
  list<T> sliced_list<T, S>::operator+(array_base<T, N, V> const &s) const
  {
    list<T> out(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), out.begin()));
    return out;
  }
  template <class T, class S>
  template <class Tp, class Sp>
  list<typename __combined<T, Tp>::type>
  sliced_list<T, S>::operator+(sliced_list<Tp, Sp> const &s) const
  {
    list<typename __combined<T, Tp>::type> out(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), out.begin()));
    return out;
  }
  template <class N, class T>
  list<T> operator*(N n, list<T> const &l)
  {
    return l * n;
  }
#ifdef USE_XSIMD
  template <class T, class S>
  template <class vectorizer>
  typename sliced_list<T, S>::simd_iterator sliced_list<T, S>::vbegin(vectorizer) const
  {
    return {_data->data() + slicing.lower};
  }

  template <class T, class S>
  template <class vectorizer>
  typename sliced_list<T, S>::simd_iterator sliced_list<T, S>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {_data->data() + slicing.lower + long(size() / vector_size * vector_size)};
  }

#endif

  // other operations
  template <class T, class S>
  template <class V>
  bool sliced_list<T, S>::contains(V const &v) const
  {
    return std::find(_data->begin(), _data->end(), v) != _data->end();
  }
  template <class T, class S>
  intptr_t sliced_list<T, S>::id() const
  {
    // sharing is not implemented for sliced list
    return reinterpret_cast<intptr_t>(this);
  }

  template <class T, class S>
  long sliced_list<T, S>::count(T const &x) const
  {
    return std::count(begin(), end(), x);
  }

  /// List

  // constructors
  template <class T>
  list<T>::list() : _data(utils::no_memory())
  {
  }
  template <class T>
  template <class InputIterator>
  list<T>::list(InputIterator start, InputIterator stop) : _data()
  {
    if (std::is_same<typename std::iterator_traits<InputIterator>::iterator_category,
                     std::random_access_iterator_tag>::value)
      _data->reserve(std::distance(start, stop));
    else
      _data->reserve(DEFAULT_CAPACITY);
    std::copy(start, stop, std::back_inserter(*_data));
  }
  template <class T>
  list<T>::list(empty_list const &) : _data(0)
  {
  }
  template <class T>
  list<T>::list(size_type sz) : _data(sz)
  {
  }
  template <class T>
  list<T>::list(std::initializer_list<T> l) : _data(std::move(l))
  {
  }
  template <class T>
  list<T>::list(list<T> &&other) : _data(std::move(other._data))
  {
  }
  template <class T>
  list<T>::list(list<T> const &other) : _data(other._data)
  {
  }
  template <class T>
  template <class F>
  list<T>::list(list<F> const &other) : _data(other.size())
  {
    std::copy(other.begin(), other.end(), begin());
  }
  template <class T>
  template <class Tp, class S>
  list<T>::list(sliced_list<Tp, S> const &other) : _data(other.begin(), other.end())
  {
  }

  // operators
  template <class T>
  list<T> &list<T>::operator=(list<T> &&other)
  {
    _data = std::move(other._data);
    return *this;
  }
  template <class T>
  template <class F>
  list<T> &list<T>::operator=(list<F> const &other)
  {
    _data = utils::shared_ref<container_type>{other.size()};
    std::copy(other.begin(), other.end(), begin());
    return *this;
  }
  template <class T>
  list<T> &list<T>::operator=(list<T> const &other)
  {
    _data = other._data;
    return *this;
  }
  template <class T>
  list<T> &list<T>::operator=(empty_list const &)
  {
    _data = utils::shared_ref<container_type>();
    return *this;
  }
  template <class T>
  template <class Tp, size_t N, class V>
  list<T> &list<T>::operator=(array_base<Tp, N, V> const &other)
  {
    _data = utils::shared_ref<container_type>(other.begin(), other.end());
    return *this;
  }
  template <class T>
  template <class Tp, class S>
  list<T> &list<T>::operator=(sliced_list<Tp, S> const &other)
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
  list<T> &list<T>::operator+=(sliced_list<T, S> const &other)
  {
    _data->resize(size() + other.size());
    std::copy(other.begin(), other.end(), _data->begin());
    return *this;
  }

  template <class T>
  template <class S>
  list<T> list<T>::operator+(sliced_list<T, S> const &other) const
  {
    list<T> new_list(begin(), end());
    new_list.reserve(size() + other.size());
    std::copy(other.begin(), other.end(), std::back_inserter(new_list));
    return new_list;
  }

  template <class T>
  template <size_t N, class V>
  list<T> list<T>::operator+(array_base<T, N, V> const &other) const
  {
    list<T> new_list(begin(), end());
    new_list.reserve(size() + other.size());
    std::copy(other.begin(), other.end(), std::back_inserter(new_list));
    return new_list;
  }

  // io
  template <class T>
  std::ostream &operator<<(std::ostream &os, list<T> const &v)
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
  bool list<T>::operator==(list<K> const &other) const
  {
    if (size() != other.size())
      return false;
    return std::equal(begin(), end(), other.begin());
  }
  template <class T>
  bool list<T>::operator==(empty_list const &) const
  {
    return size() == 0;
  }
  template <class T>
  template <class K>
  bool list<T>::operator!=(list<K> const &other) const
  {
    return !operator==(other);
  }
  template <class T>
  bool list<T>::operator!=(empty_list const &) const
  {
    return size() != 0;
  }

  // iterators
  template <class T>
  typename list<T>::iterator list<T>::begin()
  {
    return _data->begin();
  }
  template <class T>
  typename list<T>::const_iterator list<T>::begin() const
  {
    return _data->begin();
  }
  template <class T>
  typename list<T>::iterator list<T>::end()
  {
    return _data->end();
  }
  template <class T>
  typename list<T>::const_iterator list<T>::end() const
  {
    return _data->end();
  }
  template <class T>
  typename list<T>::reverse_iterator list<T>::rbegin()
  {
    return _data->rbegin();
  }
  template <class T>
  typename list<T>::const_reverse_iterator list<T>::rbegin() const
  {
    return _data->rbegin();
  }
  template <class T>
  typename list<T>::reverse_iterator list<T>::rend()
  {
    return _data->rend();
  }
  template <class T>
  typename list<T>::const_reverse_iterator list<T>::rend() const
  {
    return _data->rend();
  }

  // comparison
  template <class T>
  bool list<T>::operator<(list<T> const &other) const
  {
    return std::lexicographical_compare(begin(), end(), other.begin(), other.end());
  }
  template <class T>
  bool list<T>::operator>(list<T> const &other) const
  {
    return std::lexicographical_compare(other.begin(), other.end(), begin(), end());
  }
  template <class T>
  bool list<T>::operator<=(list<T> const &other) const
  {
    return !(*this > other);
  }
  template <class T>
  bool list<T>::operator>=(list<T> const &other) const
  {
    return !(*this < other);
  }

// element access
#ifdef USE_XSIMD
  template <class T>
  template <class vectorizer>
  typename list<T>::simd_iterator list<T>::vbegin(vectorizer) const
  {
    return {_data->data()};
  }

  template <class T>
  template <class vectorizer>
  typename list<T>::simd_iterator list<T>::vend(vectorizer) const
  {
    using vector_type = typename xsimd::batch<dtype>;
    static const std::size_t vector_size = vector_type::size;
    return {_data->data() + long(size() / vector_size * vector_size)};
  }

#endif
  template <class T>
  typename list<T>::reference list<T>::fast(long n)
  {
    return (*_data)[n];
  }
  template <class T>
  typename list<T>::reference list<T>::operator[](long n)
  {
    if (n < 0)
      n += size();
    assert(0 <= n && n < size());
    return fast(n);
  }
  template <class T>
  typename list<T>::const_reference list<T>::fast(long n) const
  {
    assert(n < size());
    return (*_data)[n];
  }
  template <class T>
  typename list<T>::const_reference list<T>::operator[](long n) const
  {
    if (n < 0)
      n += size();
    assert(0 <= n && n < size());
    return fast(n);
  }

  template <class T>
  template <class Sp>
  std::enable_if_t<is_slice<Sp>::value, sliced_list<T, Sp>> list<T>::operator[](Sp const &s) const
  {
    return {*this, s};
  }

  // modifiers
  template <class T>
  template <class Tp>
  void list<T>::push_back(Tp &&x)
  {
    // FIXME: clang-3.4 doesn't support emplace_back for vector of bool
    _data->push_back(std::forward<Tp>(x));
  }
  template <class T>
  template <class Tp>
  void list<T>::insert(long i, Tp &&x)
  {
    if (i == size())
      _data->emplace_back(std::forward<Tp>(x));
    else
      _data->insert(_data->begin() + i, std::forward<Tp>(x));
  }
  template <class T>
  void list<T>::reserve(size_t n)
  {
    if (n > _data->capacity())
      _data->reserve((n / 2) * 3);
  }
  template <class T>
  void list<T>::resize(size_t n)
  {
    _data->resize(n);
  }
  template <class T>
  typename list<T>::iterator list<T>::erase(size_t n)
  {
    return _data->erase(_data->begin() + n);
  }
  template <class T>
  T list<T>::pop(long x)
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
  void list<T>::clear()
  {
    _data->clear();
  }

  // TODO: have to raise a valueError
  template <class T>
  none_type list<T>::remove(T const &x)
  {
    erase(index(x));
    return {};
  }

  // Misc
  template <class T>
  long list<T>::index(T const &x) const
  {
    return std::find(begin(), end(), x) - begin();
  }

  // list interface
  template <class T>
  list<T>::operator bool() const
  {
    return !_data->empty();
  }

  template <class T>
  template <class F>
  list<typename __combined<T, F>::type> list<T>::operator+(list<F> const &s) const
  {
    list<typename __combined<T, F>::type> clone(size() + s.size());
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), clone.begin()));
    return clone;
  }

  template <class T>
  template <class F, class S>
  list<decltype(std::declval<T>() + std::declval<typename sliced_list<F, S>::value_type>())>
  list<T>::operator+(sliced_list<F, S> const &s) const
  {
    list<decltype(std::declval<T>() + std::declval<typename sliced_list<F, S>::value_type>())>
        clone(size() + len(s));
    std::copy(s.begin(), s.end(), std::copy(begin(), end(), clone.begin()));
    return clone;
  }

  template <class T>
  list<T> list<T>::operator+(empty_list const &) const
  {
    return list<T>(begin(), end());
  }

  template <class T>
  list<T> list<T>::operator*(long n) const
  {
    if (size() == 1) {
      list<T> r(size() * n);
      std::fill(r.begin(), r.end(), fast(0));
      return r;
    } else {
      list<T> r(size() * n);
      auto start = r.begin();
      while (start != r.end())
        start = std::copy(this->begin(), this->end(), start);
      return r;
    }
  }
  template <class T>
  list<T> const &list<T>::operator*=(long n)
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
  list<T> &list<T>::operator+=(F const &s)
  {
    reserve(size() + s.size());
    std::copy(s.begin(), s.end(), std::back_inserter(*this));
    return *this;
  }

  template <class T>
  long list<T>::size() const
  {
    return _data->size();
  }
  template <class T>
  template <class E>
  long list<T>::_flat_size(E const &e, utils::int_<1>) const
  {
    return std::distance(e.begin(), e.end());
  }
  template <class T>
  template <class E, size_t L>
  long list<T>::_flat_size(E const &e, utils::int_<L>) const
  {
    return std::distance(e.begin(), e.end()) * _flat_size(e[0], utils::int_<L - 1>{});
  }
  template <class T>
  long list<T>::flat_size() const
  {
    return _flat_size(*this, utils::int_<value>{});
  }

  template <class T>
  template <class V>
  bool list<T>::contains(V const &v) const
  {
    return std::find(_data->begin(), _data->end(), v) != _data->end();
  }
  template <class T>
  intptr_t list<T>::id() const
  {
    return reinterpret_cast<intptr_t>(&(*_data));
  }

  template <class T>
  long list<T>::count(T const &x) const
  {
    return std::count(begin(), end(), x);
  }

  /// Empty list
  template <class T>
  list<T> empty_list::operator+(list<T> const &s) const
  {
    return s;
  }
  template <class T, class S>
  sliced_list<T, S> empty_list::operator+(sliced_list<T, S> const &s) const
  {
    return s;
  }
  template <class T, size_t N, class V>
  static_list<T, N> empty_list::operator+(array_base<T, N, V> const &s) const
  {
    return s.template to_array<list_version>();
  }
  inline empty_list empty_list::operator+(empty_list const &) const
  {
    return empty_list();
  }
  template <class F>
  std::enable_if_t<!is_numexpr_arg<F>::value, list<typename F::value_type>>
  empty_list::operator+(F s) const
  {
    return {s.begin(), s.end()};
  }
  inline empty_list::operator bool() const
  {
    return false;
  }

  template <class T>
  empty_list::operator list<T>() const
  {
    return list<T>(0);
  }

  inline constexpr long empty_list::size()
  {
    return 0;
  }

  inline std::ostream &operator<<(std::ostream &os, empty_list const &)
  {
    return os << "[]";
  }
} // namespace types

namespace utils
{

  template <class T, class From>
  void reserve(types::list<T> &l, From const &f, typename From::const_iterator *)
  {
    l.reserve(builtins::len(f));
  }
} // namespace utils
PYTHONIC_NS_END

/* overload std::get */
namespace std
{
  template <size_t I, class T>
  typename pythonic::types::list<T>::reference get(pythonic::types::list<T> &t)
  {
    return t[I];
  }

  template <size_t I, class T>
  typename pythonic::types::list<T>::const_reference get(pythonic::types::list<T> const &t)
  {
    return t[I];
  }

  template <size_t I, class T>
  typename pythonic::types::list<T>::value_type get(pythonic::types::list<T> &&t)
  {
    return std::move(t)[I];
  }

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_list<T, S>::reference get(pythonic::types::sliced_list<T, S> &t)
  {
    return t[I];
  }

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_list<T, S>::const_reference
  get(pythonic::types::sliced_list<T, S> const &t)
  {
    return t[I];
  }

  template <size_t I, class T, class S>
  typename pythonic::types::sliced_list<T, S>::value_type
  get(pythonic::types::sliced_list<T, S> &&t)
  {
    return std::move(t)[I];
  }
} // namespace std

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

#ifdef Py_LIMITED_API
#define PyList_SET_ITEM PyList_SetItem
#endif

inline PyObject *to_python<typename std::vector<bool>::reference>::convert(
    typename std::vector<bool>::reference const &v)
{
  return ::to_python((bool)v);
}

inline PyObject *
to_python<std::conditional_t<std::is_same<bool, typename std::vector<bool>::const_reference>::value,
                             phantom_type, typename std::vector<bool>::const_reference>>::
    convert(typename std::vector<bool>::const_reference const &v)
{
  return ::to_python((bool)v);
}

template <class T>
PyObject *to_python<types::list<T>>::convert(types::list<T> const &v)
{
  Py_ssize_t n = v.size();
  PyObject *ret = PyList_New(n);
  for (Py_ssize_t i = 0; i < n; i++)
    PyList_SET_ITEM(ret, i, ::to_python(v[i]));
  return ret;
}
template <class T, class S>
PyObject *to_python<types::sliced_list<T, S>>::convert(types::sliced_list<T, S> const &v)
{
  Py_ssize_t n = v.size();
  PyObject *ret = PyList_New(n);
  for (Py_ssize_t i = 0; i < n; i++)
    PyList_SET_ITEM(ret, i, ::to_python(v[i]));
  return ret;
}

inline PyObject *to_python<types::empty_list>::convert(types::empty_list const &)
{
  return PyList_New(0);
}

template <class T>
bool from_python<types::list<T>>::is_convertible(PyObject *obj)
{
  if (!PyList_Check(obj))
    return false;
  if (PyObject_Not(obj))
    return true;
#ifdef Py_LIMITED_API
  PyObject *Item = PySequence_GetItem(obj, 0);
  bool result = ::is_convertible<T>(Item);
  Py_DECREF(Item);
  return result;
#else
  return ::is_convertible<T>(PySequence_Fast_GET_ITEM(obj, 0));
#endif
}

template <class T>
types::list<T> from_python<types::list<T>>::convert(PyObject *obj)
{
#ifdef Py_LIMITED_API
  Py_ssize_t l = PySequence_Size(obj);
#else
  Py_ssize_t l = PySequence_Fast_GET_SIZE(obj);
#endif
  types::list<T> v(l);

#ifdef Py_LIMITED_API
  for (Py_ssize_t i = 0; i < l; ++i) {
    PyObject *item = PySequence_GetItem(obj, i);
    v.fast(i) = ::from_python<T>(item);
    Py_DECREF(item);
  }
#else
  PyObject **core = PySequence_Fast_ITEMS(obj);
  std::transform(core, core + l, v.begin(), [](PyObject *o) { return ::from_python<T>(o); });
#endif

  return v;
}

#ifdef Py_LIMITED_API
#undef PyList_SET_ITEM
#endif

PYTHONIC_NS_END

#endif

#endif
