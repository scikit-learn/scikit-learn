#ifndef PYTHONIC_NUMPY_FLOAT_HPP
#define PYTHONIC_NUMPY_FLOAT_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/numpy/float_.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME float_
#define NUMPY_NARY_FUNC_SYM details::float64
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::float_>::convert(numpy::functor::float_ const &c)
{
  return (PyObject *)&PyDoubleArrType_Type;
}

inline bool from_python<numpy::functor::float_>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyDoubleArrType_Type;
}

inline numpy::functor::float_ from_python<numpy::functor::float_>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
