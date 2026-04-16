#ifndef PYTHONIC_NUMPY_INT8_HPP
#define PYTHONIC_NUMPY_INT8_HPP

#include "pythonic/include/numpy/int8.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline int8_t int8()
    {
      return int8_t();
    }

    template <class V>
    int8_t int8(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME int8
#define NUMPY_NARY_FUNC_SYM details::int8
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::int8>::convert(numpy::functor::int8 const &c)
{
  return (PyObject *)&PyInt8ArrType_Type;
}

inline bool from_python<numpy::functor::int8>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyInt8ArrType_Type;
}

inline numpy::functor::int8 from_python<numpy::functor::int8>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
