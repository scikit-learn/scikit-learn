#ifndef PYTHONIC_TYPES_FLOAT_HPP
#define PYTHONIC_TYPES_FLOAT_HPP

#include "pythonic/include/types/float.hpp"

#include "pythonic/types/attr.hpp"
#include <cstddef>

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"
#include <iostream>

PYTHONIC_NS_BEGIN

inline PyObject *to_python<long double>::convert(long double d)
{
  return PyArray_Scalar(&d, PyArray_DescrFromType(NPY_LONGDOUBLE), nullptr);
}

inline PyObject *to_python<double>::convert(double d)
{
  return PyFloat_FromDouble(d);
}
inline PyObject *to_python<float>::convert(float d)
{
  return PyArray_Scalar(&d, PyArray_DescrFromType(NPY_FLOAT), nullptr);
}

inline bool from_python<long double>::is_convertible(PyObject *obj)
{
  return PyArray_IsScalar(obj, LongDouble);
}

inline long double from_python<long double>::convert(PyObject *obj)
{
#ifdef Py_LIMITED_API
  npy_longdouble val;
  PyArray_ScalarAsCtype(obj, &val);
  return val;
#else
  return PyArrayScalar_VAL(obj, LongDouble);
#endif
}

inline bool from_python<double>::is_convertible(PyObject *obj)
{
  return PyFloat_Check(obj);
}
inline double from_python<double>::convert(PyObject *obj)
{
  return PyFloat_AsDouble(obj);
}

inline bool from_python<float>::is_convertible(PyObject *obj)
{
  return PyArray_IsScalar(obj, Float);
}
inline float from_python<float>::convert(PyObject *obj)
{
#ifdef Py_LIMITED_API
  npy_float val;
  PyArray_ScalarAsCtype(obj, &val);
  return val;
#else
  return PyArrayScalar_VAL(obj, Float);
#endif
}
PYTHONIC_NS_END

#endif

#endif
