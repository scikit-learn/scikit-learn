#ifndef PYTHONIC_TYPES_SET_HPP
#define PYTHONIC_TYPES_SET_HPP

#include "pythonic/include/types/set.hpp"

#include "pythonic/types/assignable.hpp"
#include "pythonic/types/empty_iterator.hpp"
#include "pythonic/types/list.hpp"

#include "pythonic/utils/iterator.hpp"
#include "pythonic/utils/reserve.hpp"
#include "pythonic/utils/shared_ref.hpp"

#include "pythonic/builtins/in.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <set>
#include <utility>

PYTHONIC_NS_BEGIN

namespace types
{

  /// set implementation
  // constructors
  template <class T>
  set<T>::set() : data(utils::no_memory())
  {
  }

  template <class T>
  template <class InputIterator>
  set<T>::set(InputIterator start, InputIterator stop) : data()
  {
    std::copy(start, stop, std::back_inserter(*this));
  }

  template <class T>
  set<T>::set(empty_set const &) : data()
  {
  }

  template <class T>
  set<T>::set(std::initializer_list<value_type> l) : data(std::move(l))
  {
  }

  template <class T>
  set<T>::set(set<T> const &other) : data(other.data)
  {
  }

  template <class T>
  template <class F>
  set<T>::set(set<F> const &other) : data()
  {
    std::copy(other.begin(), other.end(), std::inserter(*data, data->begin()));
  }

  // iterators
  template <class T>
  typename set<T>::iterator set<T>::begin()
  {
    return data->begin();
  }

  template <class T>
  typename set<T>::const_iterator set<T>::begin() const
  {
    return data->begin();
  }

  template <class T>
  typename set<T>::iterator set<T>::end()
  {
    return data->end();
  }

  template <class T>
  typename set<T>::const_iterator set<T>::end() const
  {
    return data->end();
  }

  // modifiers
  template <class T>
  T set<T>::pop()
  {
    if (size() <= 0)
      throw std::out_of_range("Trying to pop() an empty set.");

    T tmp = *begin();
    data->erase(begin());
    return tmp;
  }

  template <class T>
  void set<T>::add(const T &x)
  {
    data->insert(x);
  }

  template <class T>
  void set<T>::push_back(const T &x)
  {
    data->insert(x);
  }

  template <class T>
  void set<T>::clear()
  {
    data->clear();
  }

  template <class T>
  template <class U>
  void set<T>::discard(U const &elem)
  {
    // Remove element elem from the set if it is present.
    data->erase(elem);
  }

  template <class T>
  template <class U>
  void set<T>::remove(U const &elem)
  {
    // Remove element elem from the set. Raises KeyError if elem is !
    // contained in the set.
    if (!data->erase(elem))
      throw std::runtime_error("set.delete() : couldn't delete element ! in the set.");
  }

  // set interface
  template <class T>
  set<T>::operator bool() const
  {
    return !data->empty();
  }

  template <class T>
  long set<T>::size() const
  {
    return data->size();
  }

  // Misc

  template <class T>
  set<T> set<T>::copy() const
  {
    return set<T>(begin(), end());
  }

  template <class T>
  template <class U>
  bool set<T>::isdisjoint(U const &other) const
  {
    // Return true if the this has no elements in common with other.
    for (const_iterator it = begin(); it != end(); ++it) {
      if (in(other, *it))
        return false;
    }
    return true;
  }

  template <class T>
  template <class U>
  bool set<T>::issubset(U const &other) const
  {
    // Test whether every element in the set is in other.
    for (const_iterator it = begin(); it != end(); ++it) {
      if (!in(other, *it))
        return false;
    }
    return true;
  }

  template <class T>
  template <class U>
  bool set<T>::issuperset(U const &other) const
  {
    // Test whether every element in other is in the set.
    return other.issubset(*this);
  }

  template <class T>
  set<T> set<T>::union_() const
  {
    return set<T>(begin(), end());
  }

  template <class T>
  template <typename U, typename... Types>
  typename __combined<set<T>, U, Types...>::type set<T>::union_(U &&other, Types &&...others) const
  {
    typename __combined<set<T>, U, Types...>::type tmp = union_(std::forward<Types...>(others)...);
    tmp.data->insert(other.begin(), other.end());
    return tmp;
  }

  template <class T>
  template <typename... Types>
  none_type set<T>::update(Types &&...others)
  {
    *this = union_(std::forward<Types>(others)...);
    return {};
  }

  template <class T>
  set<T> set<T>::intersection() const
  {
    return set<T>(begin(), end());
  }

  template <class T>
  template <typename U, typename... Types>
  typename __combined<set<T>, U, Types...>::type set<T>::intersection(U const &other,
                                                                      Types const &...others) const
  {
    // Return a new set with elements common to the set && all others.
    typename __combined<set<T>, U, Types...>::type tmp = intersection(others...);
    for (auto it = begin(); it != end(); ++it) {
      if (!in(other, *it))
        tmp.discard(*it); // faster than remove() but ! direct interaction with data
    }
    return tmp;
  }

  template <class T>
  template <typename... Types>
  void set<T>::intersection_update(Types const &...others)
  {
    *this = intersection(others...);
  }

  template <class T>
  set<T> set<T>::difference() const
  {
    return set<T>(begin(), end());
  }

  template <class T>
  template <typename U, typename... Types>
  set<T> set<T>::difference(U const &other, Types const &...others) const
  {
    // Return a new set with elements in the set that are ! in the others.
    set<T> tmp = difference(others...);
    /*
                   for(iterator it=tmp.begin(); it!=tmp.end();++it){
                   if(other.get_data().find(*it)!=other.end())
                   tmp.discard(*it);
                   }
                   */ // This algo will do several times the same find(), because
       // std::set::erase() calls find. Lame!
    for (typename U::const_iterator it = other.begin(); it != other.end(); ++it) {
      tmp.discard(*it);
    }
    return tmp;
  }

  template <class T>
  template <class V>
  bool set<T>::contains(V const &v) const
  {
    return data->find(v) != data->end();
  }

  template <class T>
  template <typename... Types>
  void set<T>::difference_update(Types const &...others)
  {
    *this = difference(others...);
  }

  template <class T>
  template <typename U>
  set<typename __combined<T, U>::type> set<T>::symmetric_difference(set<U> const &other) const
  {
    // Return a new set with elements in either the set || other but ! both.
    // return ((*this-other) | (other-*this));

    // We must use fcts && ! operators because fcts have to handle any
    // iterable objects && operators only sets (cf python ref)
    return (this->difference(other)).union_(other.difference(*this));
  }

  template <class T>
  template <typename U>
  typename __combined<U, set<T>>::type set<T>::symmetric_difference(U const &other) const
  {
    // Return a new set with elements in either the set || other but ! both.
    set<typename std::iterator_traits<typename U::iterator>::value_type> tmp(other.begin(),
                                                                             other.end());

    // We must use fcts && ! operators because fcts have to handle any
    // iterable objects && operators only sets (cf python ref)
    return (this->difference(other)).union_(tmp.difference(*this));
  }

  template <class T>
  template <typename U>
  void set<T>::symmetric_difference_update(U const &other)
  {
    *this = symmetric_difference(other);
  }

  // Operators
  template <class T>
  template <class U>
  bool set<T>::operator==(set<U> const &other) const
  {
    return *data == *other.data;
  }

  template <class T>
  template <class U>
  bool set<T>::operator<=(set<U> const &other) const
  {
    // Every element in *this is in other
    return issubset(other);
  }

  template <class T>
  template <class U>
  bool set<T>::operator<(set<U> const &other) const
  {
    // Every element in this is in other && this != other
    return (*this <= other) && (this->size() != other.size());
  }

  template <class T>
  template <class U>
  bool set<T>::operator>=(set<U> const &other) const
  {
    // Every element in other is in set
    return other <= *this;
  }

  template <class T>
  template <class U>
  bool set<T>::operator>(set<U> const &other) const
  {
    // Every element in other is in set && this != other
    return other < *this;
  }

  template <class T>
  template <class U>
  set<typename __combined<T, U>::type> set<T>::operator|(set<U> const &other) const
  {
    return union_(other);
  }

  template <class T>
  template <class U>
  void set<T>::operator|=(set<U> const &other)
  {
    update(other);
  }

  template <class T>
  template <class U>
  set<typename __combined<U, T>::type> set<T>::operator&(set<U> const &other) const
  {
    return intersection(other);
  }

  template <class T>
  template <class U>
  void set<T>::operator&=(set<U> const &other)
  {
    return intersection_update(other);
  }

  template <class T>
  template <class U>
  set<T> set<T>::operator-(set<U> const &other) const
  {
    return difference(other);
  }

  template <class T>
  template <class U>
  void set<T>::operator-=(set<U> const &other)
  {
    return difference_update(other);
  }

  template <class T>
  template <class U>
  set<typename __combined<U, T>::type> set<T>::operator^(set<U> const &other) const
  {
    return symmetric_difference(other);
  }

  template <class T>
  template <class U>
  void set<T>::operator^=(set<U> const &other)
  {
    return symmetric_difference_update(other);
  }

  template <class T>
  intptr_t set<T>::id() const
  {
    return reinterpret_cast<intptr_t>(&(*data));
  }

  template <class T>
  std::ostream &operator<<(std::ostream &os, set<T> const &v)
  {
    if (v.size() == 0) {
      return os << "set()";
    }
    os << "{";
    const char *commaSeparator = "";
    for (const auto &e : v) {
      os << commaSeparator << e;
      commaSeparator = ", ";
    }
    return os << "}";
  }

  /// empty_set implementation

  inline empty_set empty_set::operator|(empty_set const &)
  {
    return empty_set();
  }

  template <class T>
  set<T> empty_set::operator|(set<T> const &s)
  {
    return s;
  }

  template <class U>
  U empty_set::operator&(U const &s)
  {
    return {};
  }

  template <class U>
  U empty_set::operator-(U const &s)
  {
    return {};
  }

  inline empty_set empty_set::operator^(empty_set const &)
  {
    return empty_set();
  }

  template <class T>
  set<T> empty_set::operator^(set<T> const &s)
  {
    return s;
  }

  template <class... Types>
  none_type empty_set::update(Types &&...)
  {
    return {};
  }

  inline empty_set::operator bool()
  {
    return false;
  }

  inline empty_set::iterator empty_set::begin() const
  {
    return empty_iterator();
  }

  inline empty_set::iterator empty_set::end() const
  {
    return empty_iterator();
  }

  template <class V>
  bool empty_set::contains(V const &) const
  {
    return false;
  }
} // namespace types
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

template <typename T>
PyObject *to_python<types::set<T>>::convert(types::set<T> const &v)
{
  PyObject *obj = PySet_New(nullptr);
  for (auto const &e : v)
    PySet_Add(obj, ::to_python(e));
  return obj;
}

inline PyObject *to_python<types::empty_set>::convert(types::empty_set)
{
  return PySet_New(nullptr);
}

template <class T>
bool from_python<types::set<T>>::is_convertible(PyObject *obj)
{
  if (PySet_Check(obj)) {
    PyObject *iterator = PyObject_GetIter(obj);
    if (PyObject *item = PyIter_Next(iterator)) {
      bool res = ::is_convertible<T>(item);
      Py_DECREF(item);
      Py_DECREF(iterator);
      return res;
    } else {
      Py_DECREF(iterator);
      return true;
    }
  }
  return false;
}
template <class T>
types::set<T> from_python<types::set<T>>::convert(PyObject *obj)
{
  types::set<T> v = types::empty_set();
  // may be useful to reserve more space ?
  PyObject *iterator = PyObject_GetIter(obj);
  while (PyObject *item = PyIter_Next(iterator)) {
    v.add(::from_python<T>(item));
    Py_DECREF(item);
  }
  Py_DECREF(iterator);
  return v;
}
PYTHONIC_NS_END

#endif

#endif
