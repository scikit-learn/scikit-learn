#ifndef PYTHONIC_TYPES_BOOL_HPP
#define PYTHONIC_TYPES_BOOL_HPP

#include "pythonic/include/types/bool.hpp"

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayobject.h"

PYTHONIC_NS_BEGIN
inline PyObject *to_python<bool>::convert(bool b)
{
  if (b)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

inline bool from_python<bool>::is_convertible(PyObject *obj)
{
  return obj == Py_True || obj == Py_False || PyObject_TypeCheck(obj, &PyBoolArrType_Type);
}
inline bool from_python<bool>::convert(PyObject *obj)
{
  if (obj == Py_True)
    return true;
  else if (obj == Py_False)
    return false;
  else
    return PyObject_IsTrue(obj);
}

PYTHONIC_NS_END

#endif

#endif
