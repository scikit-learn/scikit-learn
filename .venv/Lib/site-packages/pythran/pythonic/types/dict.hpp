#ifndef PYTHONIC_TYPES_DICT_HPP
#define PYTHONIC_TYPES_DICT_HPP

#include "pythonic/include/types/dict.hpp"

#include "pythonic/builtins/KeyError.hpp"
#include "pythonic/builtins/None.hpp"
#include "pythonic/types/empty_iterator.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/iterator.hpp"
#include "pythonic/utils/reserve.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <utility>

PYTHONIC_NS_BEGIN

namespace types
{
  /// item implementation
  template <class I>
  item_iterator_adaptator<I>::item_iterator_adaptator(I const &i) : base(i)
  {
  }
  template <class I>
  typename item_iterator_adaptator<I>::value_type item_iterator_adaptator<I>::operator*() const
  {
    auto &&tmp = *base;
    ;
    return pythonic::types::make_tuple(tmp.first, tmp.second);
  }

  /// key_iterator_adaptator implementation
  template <class I>
  key_iterator_adaptator<I>::key_iterator_adaptator(I const &i) : base(i)
  {
  }

  template <class I>
  typename key_iterator_adaptator<I>::value_type key_iterator_adaptator<I>::operator*() const
  {
    return base->first;
  }

  /// value_iterator_adaptator implementation
  template <class I>
  value_iterator_adaptator<I>::value_iterator_adaptator(I const &i) : base(i)
  {
  }

  template <class I>
  typename value_iterator_adaptator<I>::value_type value_iterator_adaptator<I>::operator*() const
  {
    return base->second;
  }

  template <class D>
  dict_items<D>::dict_items()
  {
  }

  template <class D>
  dict_items<D>::dict_items(D const &d) : data(d)
  {
  }

  template <class D>
  typename dict_items<D>::iterator dict_items<D>::begin() const
  {
    return data.item_begin();
  }

  template <class D>
  typename dict_items<D>::iterator dict_items<D>::end() const
  {
    return data.item_end();
  }

  template <class D>
  long dict_items<D>::size() const
  {
    return data.size();
  }

  template <class D>
  dict_keys<D>::dict_keys()
  {
  }

  template <class D>
  dict_keys<D>::dict_keys(D const &d) : data(d)
  {
  }

  template <class D>
  typename dict_keys<D>::iterator dict_keys<D>::begin() const
  {
    return data.key_begin();
  }

  template <class D>
  typename dict_keys<D>::iterator dict_keys<D>::end() const
  {
    return data.key_end();
  }

  template <class D>
  long dict_keys<D>::size() const
  {
    return data.size();
  }

  template <class D>
  dict_values<D>::dict_values()
  {
  }

  template <class D>
  dict_values<D>::dict_values(D const &d) : data(d)
  {
  }

  template <class D>
  typename dict_values<D>::iterator dict_values<D>::begin() const
  {
    return data.value_begin();
  }

  template <class D>
  typename dict_values<D>::iterator dict_values<D>::end() const
  {
    return data.value_end();
  }

  template <class D>
  long dict_values<D>::size() const
  {
    return data.size();
  }

  template <class K, class V>
  dict<K, V>::dict() : data(utils::no_memory())
  {
  }

  template <class K, class V>
  dict<K, V>::dict(empty_dict const &) : data(DEFAULT_DICT_CAPACITY)
  {
  }

  template <class K, class V>
  dict<K, V>::dict(std::initializer_list<value_type> l) : data(l.begin(), l.end())
  {
  }

  template <class K, class V>
  dict<K, V>::dict(dict<K, V> const &other) : data(other.data)
  {
  }

  template <class K, class V>
  template <class Kp, class Vp>
  dict<K, V>::dict(dict<Kp, Vp> const &other) : data(other.data->begin(), other.data->end())
  {
  }

  template <class K, class V>
  template <class B, class E>
  dict<K, V>::dict(B begin, E end) : data(begin, end)
  {
  }

  // iterators
  template <class K, class V>
  typename dict<K, V>::iterator dict<K, V>::begin()
  {
    return typename dict<K, V>::iterator(data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::const_iterator dict<K, V>::begin() const
  {
    return key_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::iterator dict<K, V>::end()
  {
    return typename dict<K, V>::iterator(data->end());
  }

  template <class K, class V>
  typename dict<K, V>::const_iterator dict<K, V>::end() const
  {
    return key_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(data->end());
  }

  template <class K, class V>
  typename dict<K, V>::item_iterator dict<K, V>::item_begin()
  {
    return item_iterator_adaptator<typename dict<K, V>::container_type::iterator>(data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::item_const_iterator dict<K, V>::item_begin() const
  {
    return item_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::item_iterator dict<K, V>::item_end()
  {
    return item_iterator_adaptator<typename dict<K, V>::container_type::iterator>(data->end());
  }

  template <class K, class V>
  typename dict<K, V>::item_const_iterator dict<K, V>::item_end() const
  {
    return item_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->end());
  }

  template <class K, class V>
  typename dict<K, V>::key_iterator dict<K, V>::key_begin()
  {
    return key_iterator_adaptator<typename dict<K, V>::container_type::iterator>(data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::key_const_iterator dict<K, V>::key_begin() const
  {
    return key_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::key_iterator dict<K, V>::key_end()
  {
    return key_iterator_adaptator<typename dict<K, V>::container_type::iterator>(data->end());
  }

  template <class K, class V>
  typename dict<K, V>::key_const_iterator dict<K, V>::key_end() const
  {
    return key_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(data->end());
  }

  template <class K, class V>
  typename dict<K, V>::value_iterator dict<K, V>::value_begin()
  {
    return value_iterator_adaptator<typename dict<K, V>::container_type::iterator>(data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::value_const_iterator dict<K, V>::value_begin() const
  {
    return value_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->begin());
  }

  template <class K, class V>
  typename dict<K, V>::value_iterator dict<K, V>::value_end()
  {
    return value_iterator_adaptator<typename dict<K, V>::container_type::iterator>(data->end());
  }

  template <class K, class V>
  typename dict<K, V>::value_const_iterator dict<K, V>::value_end() const
  {
    return value_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->end());
  }

  // dict interface
  template <class K, class V>
  dict<K, V>::operator bool() const
  {
    return !data->empty();
  }

  template <class K, class V>
  V &dict<K, V>::operator[](K const &key) &
  {
    return fast(key);
  }

  template <class K, class V>
  V &dict<K, V>::operator[](K const &key) const &
  {
    return fast(key);
  }

  template <class K, class V>
  V &dict<K, V>::fast(K const &key) &
  {
    return (*data)[key];
  }

  template <class K, class V>
  V &dict<K, V>::fast(K const &key) const &
  {
    auto Where = data->find(key);
    if (Where == data->end())
      throw types::KeyError(key);
    return (*data)[key];
  }

  template <class K, class V>
  typename dict<K, V>::item_const_iterator dict<K, V>::find(K const &key) const
  {
    return item_iterator_adaptator<typename dict<K, V>::container_type::const_iterator>(
        data->find(key));
  }

  template <class K, class V>
  void dict<K, V>::clear()
  {
    return data->clear();
  }

  template <class K, class V>
  dict<K, V> dict<K, V>::copy() const
  {
    return dict<K, V>(this->data->begin(), this->data->end());
  }

  template <class K, class V>
  template <class W>
  typename __combined<V, W>::type dict<K, V>::get(K const &key, W d) const
  {
    auto ivalue = data->find(key);
    if (ivalue != data->end())
      return ivalue->second;
    else
      return d;
  }

  template <class K, class V>
  none<V> dict<K, V>::get(K const &key) const
  {
    auto ivalue = data->find(key);
    if (ivalue != data->end())
      return ivalue->second;
    else
      return builtins::None;
  }

  template <class K, class V>
  template <class W>
  V &dict<K, V>::setdefault(K const &key, W d)
  {
    auto ivalue = data->find(key);
    if (ivalue != data->end())
      return ivalue->second;
    else
      return (*data)[key] = d;
  }

  template <class K, class V>
  none<V> &dict<K, V>::setdefault(K const &key)
  {
    auto ivalue = data->find(key);
    if (ivalue != data->end())
      return ivalue->second;
    else
      return (*data)[key] = builtins::None;
  }

  template <class K, class V>
  template <class K0, class W0>
  void dict<K, V>::update(dict<K0, W0> const &d)
  {
    for (auto kv : *d.data)
      (*data)[kv.first] = kv.second;
  }

  template <class K, class V>
  template <class Iterable>
  void dict<K, V>::update(Iterable const &d)
  {
    for (auto kv : d)
      (*data)[std::get<0>(kv)] = std::get<1>(kv);
  }

  template <class K, class V>
  template <class W>
  typename __combined<V, W>::type dict<K, V>::pop(K const &key, W d)
  {
    auto ivalue = data->find(key);
    if (ivalue != data->end()) {
      auto tmp = ivalue->second;
      data->erase(ivalue);
      return tmp;
    } else
      return d;
  }

  template <class K, class V>
  V dict<K, V>::pop(K const &key)
  {
    auto ivalue = data->find(key);
    if (ivalue != data->end()) {
      auto tmp = ivalue->second;
      data->erase(ivalue);
      return tmp;
    } else
      throw types::KeyError(key);
  }

  template <class K, class V>
  make_tuple_t<K, V> dict<K, V>::popitem()
  {
    auto b = data->begin();
    if (b == data->end())
      throw types::KeyError("dictionnary is empty");
    else {
      auto r = *b;
      data->erase(b);
      return make_tuple_t<K, V>{r.first, r.second};
    }
  }

  template <class K, class V>
  long dict<K, V>::size() const
  {
    return data->size();
  }

  template <class K, class V>
  dict_items<dict<K, V>> dict<K, V>::items() const
  {
    return dict_items<dict<K, V>>(*this);
  }

  template <class K, class V>
  dict_keys<dict<K, V>> dict<K, V>::keys() const
  {
    return dict_keys<dict<K, V>>(*this);
  }

  template <class K, class V>
  dict_values<dict<K, V>> dict<K, V>::values() const
  {
    return dict_values<dict<K, V>>(*this);
  }

  // id interface
  template <class K, class V>
  intptr_t dict<K, V>::id() const
  {
    return reinterpret_cast<intptr_t>(&(*data));
  }

  template <class K, class V>
  template <class T>
  bool dict<K, V>::contains(T const &key) const
  {
    return data->find(key) != data->end();
  }

  template <class K, class V>
  dict<K, V> empty_dict::operator+(dict<K, V> const &s)
  {
    return s;
  }

  inline empty_dict empty_dict::operator+(empty_dict const &)
  {
    return empty_dict();
  }

  inline empty_dict::operator bool() const
  {
    return false;
  }

  inline typename empty_dict::iterator empty_dict::begin() const
  {
    return empty_iterator();
  }

  inline typename empty_dict::iterator empty_dict::end() const
  {
    return empty_iterator();
  }

  template <class V>
  bool empty_dict::contains(V const &) const
  {
    return false;
  }

  template <class K, class V>
  dict<K, V> operator+(dict<K, V> const &d, empty_dict)
  {
    return d;
  }
} // namespace types

inline std::ostream &operator<<(std::ostream &os, types::empty_dict const &)
{
  return os << "{}";
}

template <class K, class V>
std::ostream &operator<<(std::ostream &os, std::pair<K, V> const &p)
{
  os << p.first << ": ";
  return os << p.second;
}

template <class K, class V>
std::ostream &operator<<(std::ostream &os, types::dict<K, V> const &v)
{
  os << '{';
  auto iter = v.item_begin();
  if (iter != v.item_end()) {
    auto niter = iter;
    ++niter;
    while (niter != v.item_end()) {
      os << *iter << ", ";
      ++niter, ++iter;
    }
    os << *iter;
  }
  return os << '}';
}
PYTHONIC_NS_END

/* overload std::get */
namespace std
{
  template <size_t I, class K, class V>
  auto get(pythonic::types::dict<K, V> &d) -> decltype(d[I])
  {
    return d[I];
  }

  template <size_t I, class K, class V>
  auto get(pythonic::types::dict<K, V> const &d) -> decltype(d[I])
  {
    return d[I];
  }
} // namespace std
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <typename K, typename V>
PyObject *to_python<types::dict<K, V>>::convert(types::dict<K, V> const &v)
{
  PyObject *ret = PyDict_New();
  for (auto kv = v.item_begin(); kv != v.item_end(); ++kv) {
    PyObject *kobj = ::to_python(kv->first), *vobj = ::to_python(kv->second);
    PyDict_SetItem(ret, kobj, vobj);
    Py_DECREF(kobj);
    Py_DECREF(vobj);
  }
  return ret;
}

inline PyObject *to_python<types::empty_dict>::convert(types::empty_dict)
{
  return PyDict_New();
}

template <typename K, typename V>
bool from_python<types::dict<K, V>>::

    is_convertible(PyObject *obj)
{
  if (PyDict_Check(obj)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    if (PyDict_Next(obj, &pos, &key, &value)) {
      return ::is_convertible<K>(key) && ::is_convertible<V>(value);
    } else
      return true;
  }
  return false;
}

template <typename K, typename V>
types::dict<K, V> from_python<types::dict<K, V>>::convert(PyObject *obj)
{
  types::dict<K, V> v = types::empty_dict();
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  while (PyDict_Next(obj, &pos, &key, &value))
    v[::from_python<K>(key)] = ::from_python<V>(value);
  return v;
}
PYTHONIC_NS_END

#endif

#endif
