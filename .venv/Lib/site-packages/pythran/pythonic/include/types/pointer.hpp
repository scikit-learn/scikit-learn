#ifndef PYTHONIC_INCLUDE_TYPES_POINTER_HPP
#define PYTHONIC_INCLUDE_TYPES_POINTER_HPP

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T>
  struct pointer {
    T *data;

    using reference = T &;
    using const_reference = T const &;
    using value_type = T;

    reference operator[](long);
    value_type operator[](long) const;
    reference fast(long);
    value_type fast(long) const;
  };
} // namespace types
PYTHONIC_NS_END

namespace std
{
  template <size_t I, class T>
  typename pythonic::types::pointer<T>::reference get(pythonic::types::pointer<T> &t);

  template <size_t I, class T>
  typename pythonic::types::pointer<T>::value_type get(pythonic::types::pointer<T> const &t);

  template <size_t I, class T>
  typename pythonic::types::pointer<T>::value_type get(pythonic::types::pointer<T> &&t);

  template <size_t I, class T>
  struct tuple_element<I, pythonic::types::pointer<T>> {
    typedef typename pythonic::types::pointer<T>::value_type type;
  };
} // namespace std

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <typename T>
struct to_python<types::pointer<T>> {
  static PyObject *convert(types::pointer<T> const &v);
};

template <class T>
struct from_python<types::pointer<T>> {
  static bool is_convertible(PyObject *obj);
  static types::pointer<T> convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif

#endif
