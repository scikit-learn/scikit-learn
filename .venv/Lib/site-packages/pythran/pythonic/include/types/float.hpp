#ifndef PYTHONIC_INCLUDE_TYPES_FLOAT_HPP
#define PYTHONIC_INCLUDE_TYPES_FLOAT_HPP

#include "pythonic/include/types/attr.hpp"
#include <cstddef>

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<long double> {
  static PyObject *convert(long double d);
};
template <>
struct to_python<double> {
  static PyObject *convert(double d);
};
template <>
struct to_python<float> {
  static PyObject *convert(float d);
};
template <>
struct from_python<long double> {
  static bool is_convertible(PyObject *obj);
  static long double convert(PyObject *obj);
};
template <>
struct from_python<double> {
  static bool is_convertible(PyObject *obj);
  static double convert(PyObject *obj);
};
template <>
struct from_python<float> {
  static bool is_convertible(PyObject *obj);
  static float convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif

#endif
