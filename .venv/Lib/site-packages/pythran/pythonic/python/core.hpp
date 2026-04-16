#ifndef PYTHONIC_PYTHON_CORE_HPP
#define PYTHONIC_PYTHON_CORE_HPP

#ifdef ENABLE_PYTHON_MODULE

#include "Python.h"
// Python defines this for windows, and it's not needed in C++
#undef copysign

#include <sstream>
#include <type_traits>
#include <utility>

// Cython 0.29.x still uses the deprecated API, so we can't set this macro in
// this case! Also avoid redefining it if already set by the Pythran user.
#if !defined(CYTHON_ABI) && !defined(NPY_NO_DEPRECATED_API)
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include "numpy/arrayobject.h"

PYTHONIC_NS_BEGIN
template <class T>
struct to_python;

template <class T>
struct from_python;
PYTHONIC_NS_END

template <class T>
auto to_python(T &&value)
    -> decltype(pythonic::to_python<std::remove_cv_t<std::remove_reference_t<T>>>::convert(
        std::forward<T>(value)))
{
  return pythonic::to_python<std::remove_cv_t<std::remove_reference_t<T>>>::convert(
      std::forward<T>(value));
}
template <class T>
T from_python(PyObject *obj)
{
  return pythonic::from_python<T>::convert(obj);
}
template <class T>
bool is_convertible(PyObject *obj)
{
  return pythonic::from_python<T>::is_convertible(obj);
}

PYTHONIC_NS_BEGIN

namespace python
{

#ifdef Py_LIMITED_API

  struct ByteHolder {
    PyObject *holder;
    friend std::ostream &operator<<(std::ostream &os, ByteHolder const &byte_holder)
    {
      return os << (char *)(byte_holder);
    }
    operator char *() const
    {
      return PyBytes_AsString(holder);
    }
    ~ByteHolder()
    {
      Py_DECREF(holder);
    }
  };

#define PyString_AS_STRING(obj)                                                                    \
  [obj]() {                                                                                        \
    auto *str_obj = PyUnicode_AsEncodedString(obj, "ascii", "strict");                             \
    return ::pythonic::python::ByteHolder{str_obj};                                                \
  }()
#else
#ifndef PyString_AS_STRING
#define PyString_AS_STRING (char *)_PyUnicode_COMPACT_DATA
#endif
#endif

  inline void PyObject_TypePrettyPrinter(std::ostream &oss, PyObject *obj)
  {
    if (PyTuple_Check(obj)) {
      oss << '(';
#ifdef Py_LIMITED_API
      Py_ssize_t obj_size = PyTuple_Size(obj);
      // FIXME: should we propagate the error or something?
      assert(obj_size != -1);
#else
      Py_ssize_t obj_size = PyTuple_GET_SIZE(obj);
#endif
      for (long n = obj_size, i = 0; i < n; ++i) {
#ifdef Py_LIMITED_API
        PyObject *obj_item = PyTuple_GetItem(obj, i);
        assert(obj_item);
#else
        PyObject *obj_item = PyTuple_GET_ITEM(obj, i);
#endif
        PyObject_TypePrettyPrinter(oss, obj_item);
        if (i != n - 1)
          oss << ", ";
      }
      oss << ')';
    } else if (PyArray_Check(obj)) {
      auto *arr = (PyArrayObject *)obj;
      auto *descr = PyArray_DESCR(arr);
      auto *dtype = descr->typeobj;
      auto *repr = PyObject_GetAttrString((PyObject *)dtype, "__name__");
      oss << PyString_AS_STRING(repr);
      Py_DECREF(repr);

      oss << '[';
      for (int i = 0, n = PyArray_NDIM(arr); i < n; ++i) {
        oss << ':';
        if (i != n - 1)
          oss << ", ";
      }
      oss << ']';
      if ((PyArray_FLAGS(arr) & NPY_ARRAY_F_CONTIGUOUS) &&
          ((PyArray_FLAGS(arr) & NPY_ARRAY_C_CONTIGUOUS) == 0) && (PyArray_NDIM(arr) > 1)) {
        oss << " (with unsupported column-major layout)";
      } else if (PyArray_BASE(arr)) {
        oss << " (is a view)";
      } else {
        auto const *stride = PyArray_STRIDES(arr);
        auto const *dims = PyArray_DIMS(arr);
        long current_stride = PyArray_ITEMSIZE(arr);
        for (long i = PyArray_NDIM(arr) - 1; i >= 0; i--) {
          if (stride[i] != current_stride) {
            oss << " (is strided)";
            break;
          }
          current_stride *= dims[i];
        }
      }
    } else if (PyList_Check(obj)) {
      if (PyObject_Not(obj)) {
        oss << "empty list";
      } else {
#ifdef Py_LIMITED_API
        PyObject *obj_item = PySequence_GetItem(obj, 0);
#else
        PyObject *obj_item = PySequence_Fast_GET_ITEM(obj, 0);
#endif
        PyObject_TypePrettyPrinter(oss, obj_item);
#ifdef Py_LIMITED_API
        Py_DECREF(obj_item);
#endif
        oss << " list";
      }
    } else if (PySet_Check(obj)) {
      PyObject *iterator = PyObject_GetIter(obj);
      if (PyObject *item = PyIter_Next(iterator)) {
        PyObject_TypePrettyPrinter(oss, item);
        Py_DECREF(item);
        Py_DECREF(iterator);
        oss << " set";
      } else {
        Py_DECREF(iterator);
        oss << "empty set";
      }
    } else if (PyDict_Check(obj)) {
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      if (PyDict_Next(obj, &pos, &key, &value)) {
        PyObject_TypePrettyPrinter(oss, key);
        oss << ", ";
        PyObject_TypePrettyPrinter(oss, value);
        oss << " dict";
      } else
        oss << "empty dict";
    } else if (PyCapsule_CheckExact(obj)) {
      oss << PyCapsule_GetName(obj);
    } else {
      auto *repr = PyObject_GetAttrString((PyObject *)Py_TYPE(obj), "__name__");
      oss << PyString_AS_STRING(repr);
      Py_DECREF(repr);
    }
  }

  inline std::nullptr_t raise_invalid_argument(char const name[], char const alternatives[],
                                               PyObject *args, PyObject *kwargs)
  {
    std::ostringstream oss;
    oss << "Invalid call to pythranized function `" << name << '(';
#ifdef Py_LIMITED_API
    Py_ssize_t args_size = PyTuple_Size(args);
    // FIXME: should we propagate the error or something?
    assert(args_size != -1);
#else
    Py_ssize_t args_size = PyTuple_GET_SIZE(args);
#endif
    for (long n = args_size, i = 0; i < n; ++i) {
#ifdef Py_LIMITED_API
      PyObject *args_item = PyTuple_GetItem(args, i);
      assert(args_item);
#else
      PyObject *args_item = PyTuple_GET_ITEM(args, i);
#endif
      PyObject_TypePrettyPrinter(oss, args_item);
      if (i != n - 1 || (kwargs && PyDict_Size(kwargs)))
        oss << ", ";
    }

    if (kwargs) {
      PyObject *key, *value;
      Py_ssize_t pos = 0;

      for (int next = PyDict_Next(kwargs, &pos, &key, &value); next;) {
        PyObject *vrepr = PyObject_GetAttrString((PyObject *)Py_TYPE(value), "__name__");
        oss << PyString_AS_STRING(key) << '=' << PyString_AS_STRING(vrepr);
        Py_DECREF(vrepr);
        if ((next = PyDict_Next(kwargs, &pos, &key, &value)))
          oss << ", ";
      }
    }

    oss << ")'\nCandidates are:\n" << alternatives << "\n";

    PyErr_SetString(PyExc_TypeError, oss.str().c_str());
    return nullptr;
  }
} // namespace python

PYTHONIC_NS_END

#endif

#endif
