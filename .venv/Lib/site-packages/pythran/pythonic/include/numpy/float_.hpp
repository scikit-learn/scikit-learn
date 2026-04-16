#ifndef PYTHONIC_INCLUDE_NUMPY_FLOAT_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOAT_HPP

#include "pythonic/include/numpy/float64.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME float_
#define NUMPY_NARY_FUNC_SYM details::float64
#define NUMPY_NARY_EXTRA_METHOD using type = double;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::float_> {
  static PyObject *convert(numpy::functor::float_ const &c);
};

template <>
struct from_python<numpy::functor::float_> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::float_ convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif
#endif
